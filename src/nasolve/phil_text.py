"""Line-preserving reader for PHIL parameter files (.eff, .params).

NASolve sometimes has to change a few parameters of a phenix.refine parameter
file written elsewhere (for example by the Phenix GUI) and keep everything
else exactly as it was. ``scan`` finds every scope and assignment with its
dotted path and the character span of the lines it occupies; ``edit`` removes
or replaces whole items and leaves every other character untouched. Nothing
here knows Phenix's master parameters, and nothing is reformatted.

Handled: nested ``name {`` ... ``}`` scopes, ``name = value`` assignments,
values continued with a trailing backslash, "quoted" strings (with backslash
escapes) and \"\"\"triple-quoted\"\"\" strings, and # comments.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Mapping, Sequence


class PhilSyntaxError(ValueError):
    """The text is not a PHIL file this scanner can follow."""


@dataclass(frozen=True)
class PhilItem:
    kind: str          # "scope" or "assign"
    path: str          # dotted, e.g. "refinement.refine.strategy"
    start: int         # offset of the start of the line it begins on
    end: int           # offset just past its last line (newline included)
    value: str = ""    # assignments: the value, continuations joined


def _line_start(text: str, position: int) -> int:
    return text.rfind("\n", 0, position) + 1


def _line_end(text: str, position: int) -> int:
    found = text.find("\n", position)
    return len(text) if found < 0 else found + 1


def _line_number(text: str, position: int) -> int:
    return text.count("\n", 0, position) + 1


def _read_value(text: str, position: int) -> tuple[int, str]:
    """End offset and raw text of the value starting at ``position``."""
    index, length = position, len(text)
    quote: str | None = None
    while index < length:
        char = text[index]
        if quote == '"':
            if char == "\\":
                index += 2
                continue
            if char == '"':
                quote = None
            index += 1
            continue
        if quote == '"""':
            if text.startswith('"""', index):
                quote = None
                index += 3
            else:
                index += 1
            continue
        if text.startswith('"""', index):
            quote = '"""'
            index += 3
            continue
        if char == '"':
            quote = '"'
            index += 1
            continue
        if char == "#":
            return _line_end(text, index), text[position:index]
        if char == "\\" and text[index + 1:index + 2] in ("\n", "\r"):
            index = _line_end(text, index)
            continue
        if char == "\n":
            return index + 1, text[position:index]
        if char == "}":
            return index, text[position:index]
        index += 1
    if quote is not None:
        raise PhilSyntaxError(f"unterminated string starting near line {_line_number(text, position)}")
    return length, text[position:]


def _joined(raw: str) -> str:
    return re.sub(r"[ \t]*\\\r?\n\s*", " ", raw).strip()


def scan(text: str) -> list[PhilItem]:
    """Every scope and assignment of ``text``, in document order."""
    items: list[PhilItem] = []
    stack: list[tuple[str, int]] = []
    index, length = 0, len(text)
    while index < length:
        char = text[index]
        if char in " \t\r\n":
            index += 1
            continue
        if char == "#":
            index = _line_end(text, index)
            continue
        if char == "}":
            if not stack:
                raise PhilSyntaxError(f"unmatched '}}' on line {_line_number(text, index)}")
            name, start = stack.pop()
            path = ".".join([scope for scope, _ in stack] + [name])
            items.append(PhilItem("scope", path, start, _line_end(text, index)))
            index += 1
            continue
        if char in '{="':
            raise PhilSyntaxError(f"unexpected {char!r} on line {_line_number(text, index)}")
        end = index
        while end < length and text[end] not in ' \t\r\n{}=#"':
            end += 1
        name = text[index:end]
        after = end
        while after < length and text[after] in " \t":
            after += 1
        if after < length and text[after] == "{":
            stack.append((name, _line_start(text, index)))
            index = after + 1
            continue
        if after < length and text[after] == "=":
            value_end, raw = _read_value(text, after + 1)
            path = ".".join([scope for scope, _ in stack] + [name])
            items.append(PhilItem("assign", path, _line_start(text, index), value_end, _joined(raw)))
            index = value_end
            continue
        raise PhilSyntaxError(
            f"expected '{{' or '=' after {name!r} on line {_line_number(text, index)}"
        )
    if stack:
        raise PhilSyntaxError(f"scope {stack[-1][0]!r} is never closed")
    return sorted(items, key=lambda item: (item.start, -item.end))


def unquote(value: str) -> str | None:
    """The plain string of a value; None for ``None``."""
    value = value.strip()
    if value == "None":
        return None
    if len(value) >= 6 and value.startswith('"""') and value.endswith('"""'):
        return value[3:-3]
    if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
        return re.sub(r"\\(.)", r"\1", value[1:-1])
    return value


def starred(value: str) -> tuple[str, ...]:
    """The selected words of a choice value (``*word``)."""
    return tuple(word[1:] for word in value.split() if word.startswith("*"))


def inside(items: Sequence[PhilItem], scope: PhilItem) -> list[PhilItem]:
    """Items within ``scope`` (not the scope itself)."""
    return [
        item for item in items
        if item is not scope and scope.start <= item.start and item.end <= scope.end
    ]


def edit(text: str, remove: Sequence[PhilItem] = (),
         replace: Mapping[PhilItem, str] | None = None) -> str:
    """``text`` with items removed or replaced; everything else unchanged.

    An item inside a removed item is dropped with it. Replacements give the
    full new text of an item's lines.
    """
    changes: list[tuple[int, int, str]] = [(item.start, item.end, "") for item in remove]
    for item, new_text in (replace or {}).items():
        changes.append((item.start, item.end, new_text))
    changes.sort(key=lambda change: (change[0], -change[1]))
    kept: list[tuple[int, int, str]] = []
    for change in changes:
        if kept and change[0] < kept[-1][1]:
            if change[1] <= kept[-1][1] and kept[-1][2] == "":
                continue    # inside something removed
            raise PhilSyntaxError("overlapping edits")
        kept.append(change)
    pieces, position = [], 0
    for start, end, new_text in kept:
        pieces.append(text[position:start])
        pieces.append(new_text)
        position = end
    pieces.append(text[position:])
    return "".join(pieces)


__all__ = ["PhilItem", "PhilSyntaxError", "edit", "inside", "scan", "starred", "unquote"]
