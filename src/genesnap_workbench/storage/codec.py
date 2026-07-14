"""Fixed-registry JSON codec for immutable SYN records."""

from __future__ import annotations

import json
from dataclasses import fields, is_dataclass
from datetime import date, datetime
from decimal import Decimal
from enum import Enum
from pathlib import Path

from genesnap_workbench.domain import syn as syn_domain
from genesnap_workbench.domain import shrna as shrna_domain
from genesnap_workbench.domain import expression as expression_domain
from genesnap_workbench.domain import reporter as reporter_domain
from genesnap_workbench.template_engine.syn_exports import GeneratedArtifact


def _build_type_registry() -> dict[str, type]:
    registry = {"GeneratedArtifact": GeneratedArtifact}
    for module in (syn_domain, shrna_domain, expression_domain, reporter_domain):
        for name, value in vars(module).items():
            if isinstance(value, type) and (
                is_dataclass(value) or issubclass(value, Enum)
            ):
                registry[name] = value
    return registry


TYPE_REGISTRY = _build_type_registry()


def _encode(value):
    if isinstance(value, Enum):
        return {"__enum__": type(value).__name__, "value": value.value}
    if is_dataclass(value) and not isinstance(value, type):
        return {
            "__type__": type(value).__name__,
            "fields": {
                field.name: _encode(getattr(value, field.name))
                for field in fields(value)
            },
        }
    if isinstance(value, tuple):
        return {"__tuple__": [_encode(item) for item in value]}
    if isinstance(value, list):
        return [_encode(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _encode(item) for key, item in value.items()}
    if isinstance(value, datetime):
        return {"__datetime__": value.isoformat()}
    if isinstance(value, date):
        return {"__date__": value.isoformat()}
    if isinstance(value, Decimal):
        return {"__decimal__": str(value)}
    if isinstance(value, Path):
        return {"__path__": str(value)}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"不支持持久化的数据类型：{type(value).__name__}")


def _decode(value):
    if isinstance(value, list):
        return [_decode(item) for item in value]
    if not isinstance(value, dict):
        return value
    if "__datetime__" in value:
        return datetime.fromisoformat(value["__datetime__"])
    if "__date__" in value:
        return date.fromisoformat(value["__date__"])
    if "__decimal__" in value:
        return Decimal(value["__decimal__"])
    if "__path__" in value:
        return Path(value["__path__"])
    if "__tuple__" in value:
        return tuple(_decode(item) for item in value["__tuple__"])
    if "__enum__" in value:
        enum_type = TYPE_REGISTRY.get(value["__enum__"])
        if enum_type is None or not issubclass(enum_type, Enum):
            raise ValueError(f"未知枚举类型：{value['__enum__']}")
        return enum_type(value["value"])
    if "__type__" in value:
        record_type = TYPE_REGISTRY.get(value["__type__"])
        if record_type is None or not is_dataclass(record_type):
            raise ValueError(f"未知记录类型：{value['__type__']}")
        decoded_fields = {
            name: _decode(item) for name, item in value["fields"].items()
        }
        return record_type(**decoded_fields)
    return {key: _decode(item) for key, item in value.items()}


def dumps_record(value) -> str:
    return json.dumps(
        _encode(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def loads_record(payload: str):
    return _decode(json.loads(payload))
