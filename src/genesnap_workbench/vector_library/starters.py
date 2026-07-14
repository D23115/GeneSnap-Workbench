"""Load public starter vectors only after explicit sequence confirmation."""

from __future__ import annotations

from importlib.resources import as_file, files

from Bio import SeqIO

from .models import (
    RestrictionSite,
    SiteRetentionRule,
    ShRNAVectorProtocol,
    SYNVectorProtocol,
    VectorRecord,
)


PUC57_PUBLIC_URL = "https://www.snapgene.com/plasmids/basic_cloning_vectors/pUC57"
PLKO1_PURO_PUBLIC_URL = (
    "https://www.snapgene.com/plasmids/viral_expression_and_packaging_vectors/"
    "pLKO.1_puro"
)


class StarterVectorConfirmationRequired(ValueError):
    pass


def load_public_puc57_starter(
    linearization_site: str,
    *,
    user_confirmed: bool,
) -> tuple[VectorRecord, SYNVectorProtocol]:
    """Return the public pUC57 reference after the user confirms its exact use."""
    if not user_confirmed:
        raise StarterVectorConfirmationRequired(
            "必须先确认实际使用的 pUC57 与内置公开参考序列一致",
        )
    site_key = linearization_site.strip().casefold()
    sites = {
        "ecorv": RestrictionSite("EcoRV", "GATATC", 3),
        "smai": RestrictionSite("SmaI", "CCCGGG", 3),
    }
    try:
        site = sites[site_key]
    except KeyError as error:
        raise ValueError("内置 pUC57 starter 只支持 EcoRV 或 SmaI") from error

    resource = files("genesnap_workbench.resources.vectors").joinpath(
        "puc57_snapgene_public.gb",
    )
    with as_file(resource) as path:
        record = SeqIO.read(path, "genbank")
    vector = VectorRecord.from_sequence(
        vector_record_id="public-puc57-snapgene-2710",
        structural_display_name="pUC57（SnapGene 公开参考）",
        sequence=str(record.seq),
        topology="circular",
        local_aliases=("pUC57",),
        backbone_family="pUC19-derived cloning vector",
        public_source_urls=(PUC57_PUBLIC_URL,),
        public_equivalence_status="public_reference_requires_confirmation",
    )
    protocol_id = f"puc57-{site.name.casefold()}"
    protocol = SYNVectorProtocol(
        protocol_id=protocol_id,
        protocol_version_id=f"{protocol_id}-v1",
        display_name=f"pUC57-{site.name}（公开参考，未实验验证）",
        status="enabled",
        experimental_validation_status="unverified",
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        workflow_type="de_novo_gene_synthesis",
        insertion_mode="single_digest_homology_assembly",
        linearization_site=site,
        site_retention_rule=SiteRetentionRule.NOT_REQUIRED,
        release_site=None,
        homology_arm_length=20,
    )
    return vector, protocol


def load_public_plko1_puro_starter(
    *,
    user_confirmed: bool,
) -> tuple[VectorRecord, ShRNAVectorProtocol]:
    """加载公开 pLKO.1-puro；未确认时只返回待确认 protocol。"""
    resource = files("genesnap_workbench.resources.vectors").joinpath(
        "plko1_puro_snapgene_public.gb",
    )
    with as_file(resource) as path:
        record = SeqIO.read(path, "genbank")
    vector = VectorRecord.from_sequence(
        vector_record_id="public-plko1-puro-snapgene-7050",
        structural_display_name="pLKO.1-puro（SnapGene 公开参考）",
        sequence=str(record.seq),
        topology="circular",
        local_aliases=("pLKO.1-puro", "pLKO.1"),
        backbone_family="pLKO.1",
        public_source_urls=(PLKO1_PURO_PUBLIC_URL,),
        public_equivalence_status="public_reference_requires_confirmation",
    )
    protocol = ShRNAVectorProtocol(
        protocol_id="plko1-puro-agei-ecori-shrna",
        protocol_version_id="plko1-puro-agei-ecori-shrna-v1",
        display_name="pLKO.1-puro AgeI/EcoRI shRNA（公开参考）",
        status="enabled" if user_confirmed else "pending_confirmation",
        experimental_validation_status="unverified",
        vector_record_id=vector.vector_record_id,
        vector_checksum=vector.normalized_circular_sha256,
        workflow_type="shrna_knockdown",
        insertion_mode="agei_ecori_annealed_oligo",
        left_site=RestrictionSite("AgeI", "ACCGGT", 1),
        right_site=RestrictionSite("EcoRI", "GAATTC", 1),
        sequencing_primer_name="U6",
        default_target_count=3,
        default_clones_per_target=5,
    )
    return vector, protocol
