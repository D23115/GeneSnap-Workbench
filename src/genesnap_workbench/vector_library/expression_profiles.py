"""用户电脑上的表达载体与 protocol 本地库。"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from .expression import validate_expression_protocol
from .models import ExpressionVectorProtocol, VectorRecord


class ExpressionProfileIntegrityError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class ExpressionProtocolProfileSummary:
    profile_id: str
    display_name: str
    vector_name: str
    vector_checksum: str
    protocol_version_id: str
    experimental_validation_status: str
    distribution_scope: str = "local_only"


def _safe_id(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._")
    if not clean:
        raise ValueError("protocol version ID 不能为空")
    return clean


class LocalExpressionProtocolStore:
    def __init__(self, root: Path) -> None:
        self.root = Path(root)

    def profile_path(self, profile_id: str) -> Path:
        return self.root / f"{_safe_id(profile_id)}.json"

    def save_profile(
        self,
        vector: VectorRecord,
        protocol: ExpressionVectorProtocol,
    ) -> ExpressionProtocolProfileSummary:
        validation = validate_expression_protocol(vector, protocol)
        if not validation.is_valid:
            raise ExpressionProfileIntegrityError(
                "；".join(item.message for item in validation.errors),
            )
        profile_id = _safe_id(protocol.protocol_version_id)
        self.root.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": 1,
            "distribution_scope": "local_only",
            "profile_id": profile_id,
            "vector": asdict(vector),
            "protocol": asdict(protocol),
        }
        target = self.profile_path(profile_id)
        temporary = target.with_suffix(".json.tmp")
        temporary.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(target)
        return self._summary(vector, protocol, profile_id)

    def load_profile(
        self,
        profile_id: str,
    ) -> tuple[VectorRecord, ExpressionVectorProtocol]:
        path = self.profile_path(profile_id)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("schema_version") != 1:
                raise ExpressionProfileIntegrityError("不支持的本地 protocol 版本")
            vector_data = payload["vector"]
            protocol_data = payload["protocol"]
            vector = VectorRecord.from_sequence(
                vector_record_id=vector_data["vector_record_id"],
                structural_display_name=vector_data["structural_display_name"],
                sequence=vector_data["sequence"],
                topology=vector_data["topology"],
                local_aliases=tuple(vector_data.get("local_aliases", ())),
                backbone_family=vector_data.get("backbone_family", "unknown"),
                public_source_urls=tuple(vector_data.get("public_source_urls", ())),
                public_equivalence_status=vector_data.get(
                    "public_equivalence_status",
                    "unknown",
                ),
            )
            if (
                vector.normalized_circular_sha256
                != vector_data["normalized_circular_sha256"]
            ):
                raise ExpressionProfileIntegrityError("本地载体序列校验值已变化")
            protocol = ExpressionVectorProtocol(**protocol_data)
            validation = validate_expression_protocol(vector, protocol)
            if not validation.is_valid:
                raise ExpressionProfileIntegrityError(
                    "；".join(item.message for item in validation.errors),
                )
        except ExpressionProfileIntegrityError:
            raise
        except Exception as error:
            raise ExpressionProfileIntegrityError(
                f"无法读取本地表达 protocol：{path.name}",
            ) from error
        return vector, protocol

    def list_profiles(self) -> tuple[ExpressionProtocolProfileSummary, ...]:
        if not self.root.exists():
            return ()
        summaries = []
        for path in sorted(self.root.glob("*.json")):
            vector, protocol = self.load_profile(path.stem)
            summaries.append(self._summary(vector, protocol, path.stem))
        return tuple(summaries)

    @staticmethod
    def _summary(
        vector: VectorRecord,
        protocol: ExpressionVectorProtocol,
        profile_id: str,
    ) -> ExpressionProtocolProfileSummary:
        return ExpressionProtocolProfileSummary(
            profile_id=profile_id,
            display_name=protocol.display_name,
            vector_name=vector.structural_display_name,
            vector_checksum=vector.normalized_circular_sha256,
            protocol_version_id=protocol.protocol_version_id,
            experimental_validation_status=protocol.experimental_validation_status,
        )
