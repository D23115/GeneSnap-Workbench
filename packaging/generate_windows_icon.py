"""Generate a multi-resolution Windows icon from the approved source PNG."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

from PySide6.QtCore import QByteArray, QBuffer, QIODevice, Qt
from PySide6.QtGui import QImage


ICON_SIZES = (16, 20, 24, 32, 40, 48, 64, 128, 256)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    PROJECT_ROOT
    / "src"
    / "genesnap_workbench"
    / "resources"
    / "icons"
    / "genesnap_workbench.png"
)
DEFAULT_OUTPUT = DEFAULT_SOURCE.with_suffix(".ico")


def _encode_png(image: QImage) -> bytes:
    payload = QByteArray()
    buffer = QBuffer(payload)
    if not buffer.open(QIODevice.OpenModeFlag.WriteOnly):
        raise RuntimeError("无法创建图标编码缓冲区")
    if not image.save(buffer, "PNG"):
        raise RuntimeError("无法把图标尺寸编码为 PNG")
    buffer.close()
    return bytes(payload)


def build_windows_icon(source: Path, output: Path) -> None:
    image = QImage(str(source))
    if image.isNull():
        raise ValueError(f"无法读取图标源文件：{source}")
    if image.width() != image.height():
        raise ValueError("Windows 图标源文件必须是正方形")

    encoded_images: list[tuple[int, bytes]] = []
    for size in ICON_SIZES:
        scaled = image.scaled(
            size,
            size,
            Qt.AspectRatioMode.IgnoreAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        ).convertToFormat(QImage.Format.Format_ARGB32)
        encoded_images.append((size, _encode_png(scaled)))

    header_size = 6 + 16 * len(encoded_images)
    offset = header_size
    entries: list[bytes] = []
    payloads: list[bytes] = []
    for size, payload in encoded_images:
        encoded_size = 0 if size == 256 else size
        entries.append(
            struct.pack(
                "<BBBBHHII",
                encoded_size,
                encoded_size,
                0,
                0,
                1,
                32,
                len(payload),
                offset,
            )
        )
        payloads.append(payload)
        offset += len(payload)

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(
        struct.pack("<HHH", 0, 1, len(encoded_images))
        + b"".join(entries)
        + b"".join(payloads)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build_windows_icon(args.source, args.output)
    print(f"Windows icon: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
