# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path


project_root = Path(SPECPATH).parent
source_root = project_root / "src"
resources_root = source_root / "genesnap_workbench" / "resources"

a = Analysis(
    [str(source_root / "genesnap_workbench" / "app" / "main.py")],
    pathex=[str(source_root)],
    binaries=[],
    datas=[
        (
            str(resources_root / "vectors" / "puc57_snapgene_public.gb"),
            "genesnap_workbench/resources/vectors",
        ),
        (
            str(resources_root / "vectors" / "plko1_puro_snapgene_public.gb"),
            "genesnap_workbench/resources/vectors",
        ),
        (
            str(resources_root / "icons" / "genesnap_workbench.ico"),
            "genesnap_workbench/resources/icons",
        ),
        (
            str(resources_root / "icons" / "genesnap_workbench.png"),
            "genesnap_workbench/resources/icons",
        ),
    ],
    hiddenimports=[
        "Bio.Align",
        "Bio.Blast.NCBIWWW",
        "Bio.Blast.NCBIXML",
        "Bio.SeqIO.AbiIO",
        "Bio.SeqIO.FastaIO",
        "Bio.SeqIO.InsdcIO",
        "Bio.SeqIO.SnapGeneIO",
        "docx",
        "openpyxl",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter"],
    noarchive=False,
    optimize=1,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="GeneSnapWorkbench",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(resources_root / "icons" / "genesnap_workbench.ico"),
)
