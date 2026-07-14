"""Operational PySide6 desktop interface for GeneSnap Workbench."""

from __future__ import annotations

from dataclasses import replace
from datetime import date, datetime
from decimal import Decimal
from pathlib import Path
import re
from typing import Callable
from uuid import uuid4

from PySide6.QtCore import QDate, QThread, Qt, QUrl, Signal
from PySide6.QtGui import (
    QAction,
    QColor,
    QDesktopServices,
    QFont,
    QFontDatabase,
    QKeySequence,
    QPalette,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDateEdit,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressDialog,
    QPushButton,
    QScrollArea,
    QSplitter,
    QSpinBox,
    QSizePolicy,
    QStyle,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from genesnap_workbench.project_workflow.business_calendar import (
    ChinaBusinessCalendar,
)
from genesnap_workbench.project_workflow.syn_state import display_status_label
from genesnap_workbench.domain.syn import (
    PlasmidPrepStatus,
    SYNAssemblyAttemptResult,
    SYNColonyPCRResult,
    SYNRoute,
    SYNSequencingResult,
)
from genesnap_workbench.domain.shrna import BlastScreenStatus, ShRNACandidate
from genesnap_workbench.integrations.ncbi_transcripts import (
    NCBITranscriptClient,
    TranscriptCandidate,
)
from genesnap_workbench.integrations.ncbi_blast import NCBIBlastClient
from genesnap_workbench.integrations.shrna_online import (
    ShRNAOnlineDesigner,
    ShRNAOnlineDesignResult,
)
from genesnap_workbench.sequence_core.shrna import select_initial_candidates
from genesnap_workbench.sequence_core.shrna_candidates import generate_shrna_candidates
from genesnap_workbench.project_workflow.syn_materials import MaterialReadiness
from genesnap_workbench.project_workflow.syn_service import SYNWorkflowService
from genesnap_workbench.project_workflow.syn_state import SYNStateTransitionService
from genesnap_workbench.storage.syn_repository import StoredSYNProject
from genesnap_workbench.storage.shrna_repository import StoredShRNAProject
from genesnap_workbench.storage.expression_repository import StoredExpressionProject
from genesnap_workbench.storage.reporter_repository import StoredReporterProject
from genesnap_workbench.vector_library.comparison import read_sequence_file
from genesnap_workbench.vector_library.expression import validate_expression_protocol
from genesnap_workbench.vector_library.reporter import validate_reporter_protocol
from genesnap_workbench.vector_library.models import (
    ExpressionVectorProtocol,
    ReporterVectorProtocol,
    VectorRecord,
)
from genesnap_workbench.template_engine.workbook_templates import (
    ContactProfile,
    WorkbookTemplateInspection,
    contact_field_names,
    inspect_workbook_template,
    table_field_names,
    workbook_mapping_choices,
)

from .application import (
    GeneSnapApplicationService,
    NewExpressionProjectCommand,
    NewReporterProjectCommand,
    NewSYNProjectCommand,
    NewShRNAProjectCommand,
    PreparedSYNProject,
)


APP_STYLESHEET = """
QWidget {
    color: #1f2a2d;
    background: #f5f7f7;
    font-family: "Microsoft YaHei UI", "Segoe UI";
    font-size: 13px;
}
QMainWindow, QDialog { background: #f5f7f7; }
QLabel { background: transparent; }
QAbstractScrollArea,
QAbstractScrollArea QWidget#qt_scrollarea_viewport {
    background: #ffffff;
}
QTableView, QTableWidget {
    background: #ffffff;
    alternate-background-color: #f4f7f7;
}
QTabWidget::pane { background: #ffffff; }
QToolBar {
    background: #ffffff;
    border: 0;
    border-bottom: 1px solid #d8dede;
    spacing: 6px;
    padding: 6px 10px;
}
QToolButton, QPushButton {
    background: #ffffff;
    border: 1px solid #b8c4c3;
    border-radius: 4px;
    min-height: 28px;
    padding: 2px 10px;
}
QToolButton:hover, QPushButton:hover { background: #edf5f3; }
QToolButton:pressed, QPushButton:pressed { background: #dfecea; }
QPushButton#primaryButton {
    background: #176b5b;
    border-color: #176b5b;
    color: #ffffff;
    font-weight: 600;
}
QPushButton#primaryButton:hover { background: #125a4d; }
QLineEdit, QPlainTextEdit, QComboBox, QDateEdit {
    background: #ffffff;
    border: 1px solid #b8c4c3;
    border-radius: 3px;
    padding: 5px 7px;
    selection-background-color: #8bc4b9;
}
QLineEdit:focus, QPlainTextEdit:focus, QComboBox:focus, QDateEdit:focus {
    border-color: #176b5b;
}
QTableWidget {
    background: #ffffff;
    alternate-background-color: #f7f9f9;
    border: 1px solid #d8dede;
    gridline-color: #e4e9e9;
    selection-background-color: #cfe5df;
    selection-color: #17201f;
}
QHeaderView::section {
    background: #eef2f2;
    border: 0;
    border-right: 1px solid #d8dede;
    border-bottom: 1px solid #cbd4d4;
    padding: 7px 8px;
    font-weight: 600;
}
QTabWidget::pane { background: #ffffff; border: 1px solid #d8dede; }
QTabBar::tab {
    background: #e9eeee;
    border: 1px solid #d1d9d9;
    padding: 7px 14px;
    margin-right: 2px;
}
QTabBar::tab:selected { background: #ffffff; border-bottom-color: #ffffff; }
QLabel#appTitle { font-size: 18px; font-weight: 700; color: #153b35; }
QLabel#detailTitle { font-size: 16px; font-weight: 700; color: #153b35; }
QLabel#mutedLabel { color: #667573; }
QFrame#detailHeader { background: #ffffff; border: 1px solid #d8dede; }
"""


def _lookup_transcript_candidate(
    parent: QWidget,
    client: NCBITranscriptClient,
    *,
    accession: str,
    gene_symbol: str,
    species: str,
) -> TranscriptCandidate | None:
    try:
        if accession.strip():
            return client.fetch_accession(accession.strip())
        if not gene_symbol.strip():
            raise ValueError("请先填写基因名，或直接填写转录本号")
        candidates = client.find_candidates(gene_symbol.strip(), species)
    except (ValueError, LookupError, ConnectionError) as error:
        QMessageBox.warning(parent, "NCBI 查询失败", str(error))
        return None

    if species == "human":
        candidates = tuple(item for item in candidates if item.is_mane_select)
        if not candidates:
            QMessageBox.warning(
                parent,
                "未找到 MANE Select",
                "未找到可确认的 MANE Select 转录本。请填写明确的转录本号，或手动粘贴 CDS/FASTA。",
            )
            return None
        if len(candidates) == 1:
            return candidates[0]

    labels = [item.display_label for item in candidates]
    selected_label, accepted = QInputDialog.getItem(
        parent,
        "选择转录本",
        "请确认用于设计的转录本：",
        labels,
        0,
        False,
    )
    if not accepted:
        return None
    return candidates[labels.index(selected_label)]


def _populate_workbook_template_combo(
    combo: QComboBox,
    service: GeneSnapApplicationService | None,
    kind: str,
) -> None:
    combo.addItem("内置标准模板", None)
    if service is None:
        return
    for profile in service.list_workbook_templates(kind):
        combo.addItem(profile.display_name, profile.template_id)


def _vendor_name_from_combo(combo: QComboBox) -> str:
    if combo.currentData() is None:
        return "标准"
    name = combo.currentText().strip()
    for suffix in ("引物订购表", "测序订购表", "测序表", "送测表"):
        if name.endswith(suffix):
            name = name[: -len(suffix)].strip()
    return name or "自定义"


def _project_clone_names(snapshot, prefix: str, initial_count: int) -> tuple[str, ...]:
    names = {f"{prefix}-{clone_no}" for clone_no in range(1, initial_count + 1)}
    for submission in snapshot.sequencing_submissions:
        names.update(
            sample_name
            for sample_name in submission.sample_names
            if sample_name.startswith(f"{prefix}-")
        )
    return tuple(
        sorted(
            names,
            key=lambda value: (
                int(value.rsplit("-", 1)[-1])
                if value.rsplit("-", 1)[-1].isdigit()
                else 10**9,
                value,
            ),
        ),
    )


def _sequencing_tracking_rows(snapshot) -> tuple[tuple[str, str], ...]:
    if not snapshot.sequencing_submissions:
        return ()
    latest = snapshot.sequencing_submissions[-1]
    history = "；".join(
        (
            f"第{submission.round_no}轮 {submission.submission_kind} "
            f"{submission.status} / {len(submission.sample_names)}样本 / "
            f"送测编号:{submission.internal_submission_no or '-'} / "
            f"订单号:{submission.vendor_order_no or '-'}"
        )
        for submission in snapshot.sequencing_submissions
    )
    return (
        ("最新送测编号", latest.internal_submission_no),
        ("最新订单号", latest.vendor_order_no),
        ("送测历史", history),
    )


def ensure_chinese_font() -> None:
    """Register a Windows CJK font when Qt's platform plugin exposes none."""
    candidates = (
        Path("C:/Windows/Fonts/msyh.ttc"),
        Path("C:/Windows/Fonts/simhei.ttf"),
        Path("C:/Windows/Fonts/Deng.ttf"),
    )
    for path in candidates:
        if path.exists() and QFontDatabase.addApplicationFont(str(path)) >= 0:
            QApplication.setFont(QFont("Microsoft YaHei", 10))
            return


class CopyableTableWidget(QTableWidget):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Base, QColor("#ffffff"))
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor("#f7f9f9"))
        palette.setColor(QPalette.ColorRole.Window, QColor("#ffffff"))
        self.setPalette(palette)
        self.viewport().setAutoFillBackground(True)

    def headers(self) -> tuple[str, ...]:
        return tuple(
            self.horizontalHeaderItem(column).text()
            for column in range(self.columnCount())
        )

    def copy_selected_cells(self) -> None:
        indexes = sorted(
            self.selectedIndexes(),
            key=lambda index: (index.row(), index.column()),
        )
        if not indexes:
            return
        rows: dict[int, dict[int, str]] = {}
        for index in indexes:
            rows.setdefault(index.row(), {})[index.column()] = str(index.data() or "")
        columns = sorted({index.column() for index in indexes})
        text = "\n".join(
            "\t".join(rows[row].get(column, "") for column in columns)
            for row in sorted(rows)
        )
        QApplication.clipboard().setText(text)

    def keyPressEvent(self, event) -> None:
        if event.matches(QKeySequence.StandardKey.Copy):
            self.copy_selected_cells()
            event.accept()
            return
        super().keyPressEvent(event)


def _project_intake_layout(
    dialog: QDialog,
    *,
    preferred_width: int,
    preferred_height: int,
) -> tuple[QVBoxLayout, QVBoxLayout]:
    """Keep long project forms usable on scaled or short Windows displays."""
    screen = dialog.screen() or QApplication.primaryScreen()
    width = preferred_width
    height = preferred_height
    if screen is not None:
        available = screen.availableGeometry()
        width = min(width, max(480, available.width() - 64))
        height = min(height, max(420, available.height() - 64))

    dialog.setMinimumSize(min(480, width), min(420, height))
    dialog.resize(width, height)

    root_layout = QVBoxLayout(dialog)
    scroll_area = QScrollArea(dialog)
    scroll_area.setObjectName("projectIntakeScrollArea")
    scroll_area.setWidgetResizable(True)
    scroll_area.setFrameShape(QFrame.Shape.NoFrame)
    scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
    scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)

    content = QWidget(scroll_area)
    content_layout = QVBoxLayout(content)
    content_layout.setContentsMargins(0, 0, 4, 0)
    scroll_area.setWidget(content)
    root_layout.addWidget(scroll_area, 1)
    return root_layout, content_layout


def _configure_project_date_edits(
    received_date: QDateEdit,
    due_date: QDateEdit,
) -> None:
    for editor in (received_date, due_date):
        editor.setCalendarPopup(True)
        editor.setDisplayFormat("yyyy-MM-dd")
    due_date.setMinimumDate(received_date.date())


def _set_default_due_date(
    calendar: ChinaBusinessCalendar,
    received_date: QDateEdit,
    due_date: QDateEdit,
    workdays: int,
) -> None:
    received = received_date.date().toPython()
    calculated = calendar.add_workdays(received, workdays)
    due_date.setMinimumDate(received_date.date())
    due_date.setDate(QDate(calculated.year, calculated.month, calculated.day))


class NewSYNProjectDialog(QDialog):
    def __init__(
        self,
        calendar: ChinaBusinessCalendar,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.calendar = calendar
        self.setWindowTitle("新建 SYN 项目")
        layout, content_layout = _project_intake_layout(
            self,
            preferred_width=680,
            preferred_height=640,
        )
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self.project_id = QLineEdit()
        self.target_name = QLineEdit()
        self.site = QComboBox()
        self.site.addItems(("EcoRV", "SmaI"))
        self.received_date = QDateEdit(QDate.currentDate())
        self.due_date = QDateEdit()
        _configure_project_date_edits(self.received_date, self.due_date)
        self.received_date.dateChanged.connect(self._update_due_date)
        self.sequence = QPlainTextEdit()
        self.sequence.setPlaceholderText("DNA 或单条 FASTA")
        self.sequence.setMinimumHeight(270)
        self.vector_confirmation = QCheckBox(
            "确认实际 pUC57 与内置 SnapGene 公开参考序列一致",
        )
        form.addRow("项目号", self.project_id)
        form.addRow("目标名称", self.target_name)
        form.addRow("线性化位点", self.site)
        form.addRow("接收日期", self.received_date)
        form.addRow("标准完工日期", self.due_date)
        form.addRow("目标序列", self.sequence)
        form.addRow("载体确认", self.vector_confirmation)
        content_layout.addLayout(form)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("预计算设计")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._validate_and_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._update_due_date()

    def _update_due_date(self) -> None:
        _set_default_due_date(
            self.calendar,
            self.received_date,
            self.due_date,
            15,
        )

    def _validate_and_accept(self) -> None:
        if not self.project_id.text().strip():
            QMessageBox.warning(self, "缺少项目号", "项目号不能为空。")
            self.project_id.setFocus()
            return
        if not self.target_name.text().strip():
            QMessageBox.warning(self, "缺少目标名称", "目标名称不能为空。")
            self.target_name.setFocus()
            return
        if not self.sequence.toPlainText().strip():
            QMessageBox.warning(self, "缺少序列", "目标 DNA/FASTA 不能为空。")
            self.sequence.setFocus()
            return
        if not self.vector_confirmation.isChecked():
            QMessageBox.warning(self, "需要载体确认", "请先确认实际 pUC57 序列。")
            return
        self.accept()

    def command(self, actor: str) -> NewSYNProjectCommand:
        received = self.received_date.date()
        due = self.due_date.date()
        return NewSYNProjectCommand(
            project_id=self.project_id.text().strip(),
            target_name=self.target_name.text().strip(),
            raw_sequence=self.sequence.toPlainText(),
            input_format=(
                "fasta"
                if self.sequence.toPlainText().lstrip().startswith(">")
                else "plain"
            ),
            linearization_site=self.site.currentText(),
            received_date=date(received.year(), received.month(), received.day()),
            due_date=date(due.year(), due.month(), due.day()),
            actor=actor,
            vector_sequence_confirmed=self.vector_confirmation.isChecked(),
        )


class NewProjectDialog(QDialog):
    """Choose the workflow before opening its detailed intake form."""

    def __init__(self, projects_root: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("新建项目")
        self.resize(620, 220)
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.project_type = QComboBox()
        self.project_type.addItem("shRNA 敲低", "shrna")
        self.project_type.addItem("表达类", "expression")
        self.project_type.addItem("GL002 报告载体", "reporter")
        self.project_type.addItem("全基因合成", "syn")
        self.output_root = QLineEdit(str(projects_root))
        output_row = QWidget()
        output_layout = QHBoxLayout(output_row)
        output_layout.setContentsMargins(0, 0, 0, 0)
        output_layout.addWidget(self.output_root, 1)
        browse = QPushButton("选择文件夹")
        browse.clicked.connect(self._choose_output_root)
        output_layout.addWidget(browse)
        self.description = QLabel()
        self.description.setWordWrap(True)
        self.project_type.currentIndexChanged.connect(self._update_description)
        form.addRow("项目类型", self.project_type)
        form.addRow("项目保存位置", output_row)
        form.addRow("说明", self.description)
        layout.addLayout(form)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("继续")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._update_description()

    def _choose_output_root(self) -> None:
        chosen = QFileDialog.getExistingDirectory(
            self,
            "选择 GeneSnap 项目保存位置",
            self.output_root.text(),
        )
        if chosen:
            self.output_root.setText(chosen)

    def _update_description(self) -> None:
        descriptions = {
            "shrna": "Broad GPP 设计 target，NCBI BLAST 筛选，并生成 pLKO 引物与送测文件。",
            "expression": "全长 CDS、截短体、缺失体和点突变等表达构建。",
            "reporter": "GL002 promoter-luciferase 的 WT、删除体和突变体构建。",
            "syn": "客户序列 QC、overlapping oligo、两轮 PCR 与 pUC57 组装。",
        }
        self.description.setText(descriptions[str(self.project_type.currentData())])

    def _accept_validated(self) -> None:
        path = self.output_root.text().strip()
        if not path:
            QMessageBox.warning(self, "缺少保存位置", "请选择项目保存位置。")
            return
        try:
            Path(path).mkdir(parents=True, exist_ok=True)
        except OSError as error:
            QMessageBox.warning(self, "保存位置不可用", str(error))
            return
        self.accept()

    @property
    def workflow_type(self) -> str:
        return str(self.project_type.currentData())


class BroadTermsDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("首次使用 Broad GPP")
        self.resize(620, 300)
        layout = QVBoxLayout(self)
        notice = QLabel(
            "GeneSnap 将把当前 CDS 提交到 Broad GPP 的公开 hairpin 设计页面。"
            "Broad 当前条款说明该服务面向研究用途，商业用途可能需要另行许可。"
            "软件不会替你在网页上静默同意条款。",
        )
        notice.setWordWrap(True)
        layout.addWidget(notice)
        link = QLabel(
            '<a href="https://portals.broadinstitute.org/gpp/public/terms_and_conditions">'
            "查看 Broad GPP 服务条款</a>",
        )
        link.setOpenExternalLinks(True)
        layout.addWidget(link)
        self.confirmation = QCheckBox("我已阅读并同意按 Broad GPP 当前条款使用该在线服务")
        layout.addWidget(self.confirmation)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("同意并继续")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept_validated(self) -> None:
        if not self.confirmation.isChecked():
            QMessageBox.warning(self, "尚未确认", "请先阅读并确认 Broad GPP 服务条款。")
            return
        self.accept()


class ShRNAOnlineDesignThread(QThread):
    progress = Signal(str)
    succeeded = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        designer: ShRNAOnlineDesigner,
        *,
        cds_sequence: str,
        gene_symbol: str,
        species: str,
        target_count: int,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.designer = designer
        self.arguments = {
            "cds_sequence": cds_sequence,
            "gene_symbol": gene_symbol,
            "species": species,
            "target_count": target_count,
        }

    def run(self) -> None:
        try:
            result = self.designer.design(progress=self.progress.emit, **self.arguments)
        except Exception as error:
            self.failed.emit(str(error))
            return
        self.succeeded.emit(result)


class NewShRNAProjectDialog(QDialog):
    def __init__(
        self,
        calendar: ChinaBusinessCalendar,
        parent: QWidget | None = None,
        *,
        ncbi_client: NCBITranscriptClient | None = None,
        online_designer: ShRNAOnlineDesigner | None = None,
        service: GeneSnapApplicationService | None = None,
    ) -> None:
        super().__init__(parent)
        self.calendar = calendar
        self.ncbi_client = ncbi_client or NCBITranscriptClient()
        self.service = service
        contact_email = ""
        if service is not None:
            contact_email = service.contact_profile_store.load().email
        self.online_designer = online_designer or ShRNAOnlineDesigner(
            blast_client=NCBIBlastClient(email=contact_email),
        )
        self._online_thread: ShRNAOnlineDesignThread | None = None
        self._online_progress: QProgressDialog | None = None
        self.setWindowTitle("新建 shRNA 项目")
        layout, content_layout = _project_intake_layout(
            self,
            preferred_width=700,
            preferred_height=720,
        )
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self.project_id = QLineEdit()
        self.gene_symbol = QLineEdit()
        self.species = QComboBox()
        self.species.addItems(("human", "mouse", "rat"))
        self.transcript_accession = QLineEdit()
        transcript_row = QWidget()
        transcript_layout = QHBoxLayout(transcript_row)
        transcript_layout.setContentsMargins(0, 0, 0, 0)
        transcript_layout.addWidget(self.transcript_accession, 1)
        self.lookup_transcript_button = QPushButton("查询 CDS")
        self.lookup_transcript_button.clicked.connect(self._lookup_transcript)
        transcript_layout.addWidget(self.lookup_transcript_button)
        self.cds_sequence = QPlainTextEdit()
        self.cds_sequence.setMinimumHeight(150)
        self.targets = QPlainTextEdit()
        self.targets.setMinimumHeight(110)
        self.targets.setPlaceholderText("每行一个已经确认的 target 序列")
        self._generated_candidates_by_sequence: dict[str, ShRNACandidate] = {}
        candidate_tools = QVBoxLayout()
        candidate_buttons = QHBoxLayout()
        self.online_design_button = QPushButton("Broad + NCBI BLAST 自动设计")
        self.online_design_button.setObjectName("primaryButton")
        self.online_design_button.clicked.connect(self._start_online_design)
        self.generate_candidates_button = QPushButton("离线候选（备用）")
        self.generate_candidates_button.clicked.connect(self._generate_candidates)
        self.candidate_summary = QLabel("默认走 Broad GPP + NCBI BLAST；也可手动粘贴 target")
        self.candidate_summary.setWordWrap(True)
        candidate_buttons.addWidget(self.online_design_button)
        candidate_buttons.addWidget(self.generate_candidates_button)
        candidate_buttons.addStretch(1)
        candidate_tools.addLayout(candidate_buttons)
        candidate_tools.addWidget(self.candidate_summary)
        self.online_result_table = QTableWidget(0, 6)
        self.online_result_table.setHorizontalHeaderLabels(
            ("Target", "Broad 得分", "位置", "NCBI BLAST", "脱靶信息", "Oligo 来源"),
        )
        self.online_result_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.online_result_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.online_result_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents,
        )
        self.online_result_table.horizontalHeader().setStretchLastSection(True)
        self.online_result_table.setMinimumHeight(130)
        self.candidate_confirmation = QCheckBox(
            "已人工确认当前候选可继续设计（自动 BLAST 尚未执行）",
        )
        self.candidate_confirmation.setEnabled(False)
        self.target_count = QSpinBox()
        self.target_count.setRange(1, 3)
        self.target_count.setValue(3)
        self.clones_per_target = QSpinBox()
        self.clones_per_target.setRange(1, 96)
        self.clones_per_target.setValue(5)
        self.primer_template_combo = QComboBox()
        _populate_workbook_template_combo(
            self.primer_template_combo,
            service,
            "primer_order",
        )
        self.sequencing_template_combo = QComboBox()
        _populate_workbook_template_combo(
            self.sequencing_template_combo,
            service,
            "sequencing_order",
        )
        self.received_date = QDateEdit(QDate.currentDate())
        self.due_date = QDateEdit()
        _configure_project_date_edits(self.received_date, self.due_date)
        self.received_date.dateChanged.connect(self._update_due_date)
        self.vector_confirmation = QCheckBox(
            "已确认实际载体与所选 pLKO.1 protocol 相符",
        )
        form.addRow("项目号", self.project_id)
        form.addRow("基因", self.gene_symbol)
        form.addRow("物种", self.species)
        form.addRow("转录本号（可选）", transcript_row)
        form.addRow("CDS / FASTA", self.cds_sequence)
        form.addRow("候选生成", candidate_tools)
        form.addRow("在线筛选结果", self.online_result_table)
        form.addRow("Target", self.targets)
        form.addRow("人工确认", self.candidate_confirmation)
        form.addRow("Target 数量", self.target_count)
        form.addRow("每个 target 送测克隆数", self.clones_per_target)
        form.addRow("引物订购模板", self.primer_template_combo)
        form.addRow("送测模板", self.sequencing_template_combo)
        form.addRow("接收日期", self.received_date)
        form.addRow("标准完工日期", self.due_date)
        form.addRow("载体确认", self.vector_confirmation)
        content_layout.addLayout(form)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("生成设计并保存")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._update_due_date()

    def _update_due_date(self) -> None:
        _set_default_due_date(
            self.calendar,
            self.received_date,
            self.due_date,
            9,
        )

    def _target_lines(self) -> tuple[str, ...]:
        return tuple(
            line.strip().upper()
            for line in self.targets.toPlainText().splitlines()
            if line.strip()
        )

    def _lookup_transcript(self) -> None:
        candidate = _lookup_transcript_candidate(
            self,
            self.ncbi_client,
            accession=self.transcript_accession.text(),
            gene_symbol=self.gene_symbol.text(),
            species=self.species.currentText(),
        )
        if candidate is None:
            return
        self.transcript_accession.setText(candidate.accession)
        if candidate.gene_symbol:
            self.gene_symbol.setText(candidate.gene_symbol)
        self.cds_sequence.setPlainText(candidate.cds_sequence)

    def _generate_candidates(self) -> None:
        try:
            candidates = generate_shrna_candidates(
                self.cds_sequence.toPlainText(),
                max_candidates=60,
            )
            selection = select_initial_candidates(
                candidates,
                target_count=self.target_count.value(),
            )
        except ValueError as error:
            QMessageBox.warning(self, "无法生成候选", str(error))
            return
        self._generated_candidates_by_sequence = {
            item.target_sequence: item for item in candidates
        }
        self.targets.setPlainText(
            "\n".join(item.target_sequence for item in selection.selected),
        )
        spacing_note = "已满足 >100 bp 间距" if not selection.spacing_relaxed else "已放宽间距"
        self.candidate_summary.setText(
            f"已生成 {len(candidates)} 条候选，选出 {len(selection.selected)} 条；{spacing_note}",
        )
        self.candidate_confirmation.setEnabled(True)
        self.candidate_confirmation.setChecked(False)
        self.online_result_table.setRowCount(0)

    def _start_online_design(self) -> None:
        if not self.gene_symbol.text().strip() or not self.cds_sequence.toPlainText().strip():
            QMessageBox.warning(self, "缺少输入", "请先填写基因并查询或粘贴 CDS。")
            return
        if self.service is not None and not self.service.has_accepted_broad_terms():
            terms = BroadTermsDialog(self)
            if terms.exec() != QDialog.DialogCode.Accepted:
                return
            self.service.accept_broad_terms()
        self.online_design_button.setEnabled(False)
        self.generate_candidates_button.setEnabled(False)
        self._online_progress = QProgressDialog(
            "正在准备在线设计…",
            "",
            0,
            0,
            self,
        )
        self._online_progress.setWindowTitle("shRNA 自动设计")
        self._online_progress.setCancelButton(None)
        self._online_progress.setWindowModality(Qt.WindowModality.WindowModal)
        self._online_progress.show()
        self._online_thread = ShRNAOnlineDesignThread(
            self.online_designer,
            cds_sequence=self.cds_sequence.toPlainText(),
            gene_symbol=self.gene_symbol.text().strip(),
            species=self.species.currentText(),
            target_count=self.target_count.value(),
            parent=self,
        )
        self._online_thread.progress.connect(self._set_online_progress)
        self._online_thread.succeeded.connect(self._online_design_succeeded)
        self._online_thread.failed.connect(self._online_design_failed)
        self._online_thread.finished.connect(self._online_thread.deleteLater)
        self._online_thread.start()

    def _set_online_progress(self, message: str) -> None:
        if self._online_progress is not None:
            self._online_progress.setLabelText(message)

    def _finish_online_thread(self) -> None:
        if self._online_progress is not None:
            self._online_progress.close()
            self._online_progress = None
        self.online_design_button.setEnabled(True)
        self.generate_candidates_button.setEnabled(True)
        self._online_thread = None

    def reject(self) -> None:
        if self._online_thread is not None and self._online_thread.isRunning():
            QMessageBox.information(
                self,
                "在线设计正在运行",
                "请等待当前 Broad/NCBI 步骤结束后再关闭窗口。",
            )
            return
        super().reject()

    def closeEvent(self, event) -> None:
        if self._online_thread is not None and self._online_thread.isRunning():
            event.ignore()
            return
        super().closeEvent(event)

    def _online_design_succeeded(self, result: ShRNAOnlineDesignResult) -> None:
        self._apply_online_design_result(result)
        self._finish_online_thread()

    def _online_design_failed(self, message: str) -> None:
        self._finish_online_thread()
        QMessageBox.warning(
            self,
            "在线设计未完成",
            message + "\n\n可以重试，或使用“离线候选（备用）”后人工 BLAST。",
        )

    def _apply_online_design_result(self, result: ShRNAOnlineDesignResult) -> None:
        self._generated_candidates_by_sequence = {
            item.target_sequence: item for item in result.candidate_pool
        }
        for item in result.selected_candidates:
            self._generated_candidates_by_sequence[item.target_sequence] = item
        self.targets.setPlainText(
            "\n".join(item.target_sequence for item in result.selected_candidates),
        )
        self.online_result_table.setRowCount(len(result.selected_candidates))
        for row, item in enumerate(result.selected_candidates):
            values = (
                item.target_sequence,
                str(item.intrinsic_score),
                "" if item.start_position is None else str(item.start_position),
                item.blast_status.value,
                item.blast_note or "",
                item.oligo_source or "",
            )
            for column, value in enumerate(values):
                self.online_result_table.setItem(row, column, QTableWidgetItem(value))
        passed = sum(
            item.blast_status is BlastScreenStatus.PASS
            for item in result.selected_candidates
        )
        self.candidate_summary.setText(
            f"Broad GPP 候选；NCBI BLAST 自动通过 {passed}/{len(result.selected_candidates)}",
        )
        self.candidate_confirmation.setEnabled(result.requires_manual_confirmation)
        self.candidate_confirmation.setChecked(False)
        if result.requires_manual_confirmation:
            self.candidate_confirmation.setText(
                "已人工复核未自动通过项，并确认按当前 target/oligo 继续设计",
            )

    def _accept_validated(self) -> None:
        if not self.project_id.text().strip() or not self.gene_symbol.text().strip():
            QMessageBox.warning(self, "缺少项目信息", "项目号和基因不能为空。")
            return
        if not self.cds_sequence.toPlainText().strip():
            QMessageBox.warning(self, "缺少 CDS", "请粘贴 CDS 或 FASTA。")
            return
        lines = self._target_lines()
        if len(lines) < self.target_count.value():
            QMessageBox.warning(self, "Target 数量不足", "Target 行数少于选择的数量。")
            return
        invalid = [
            line
            for line in lines
            if not 15 <= len(line) <= 30 or set(line) - set("ACGT")
        ]
        if invalid:
            QMessageBox.warning(
                self,
                "Target 无效",
                "Target 必须是 15-30 nt 的 A/C/G/T 序列。",
            )
            return
        if self.candidate_confirmation.isEnabled() and not self.candidate_confirmation.isChecked():
            QMessageBox.warning(
                self,
                "候选尚未确认",
                "自动生成的候选尚未执行 BLAST，请先人工复核并确认后再继续。",
            )
            return
        if not self.vector_confirmation.isChecked():
            QMessageBox.warning(
                self,
                "载体尚未确认",
                "请先确认实际载体与所选 protocol 相符。",
            )
            return
        self.accept()

    def command(self, actor: str) -> NewShRNAProjectCommand:
        candidates: list[ShRNACandidate] = []
        for index, sequence in enumerate(self._target_lines(), start=1):
            generated = self._generated_candidates_by_sequence.get(sequence)
            if generated is not None:
                if generated.blast_status is BlastScreenStatus.PASS:
                    candidates.append(generated)
                else:
                    candidates.append(
                        replace(
                            generated,
                            blast_status=BlastScreenStatus.MANUALLY_ACCEPTED,
                            blast_note=(generated.blast_note or "") + "；用户已人工确认继续设计",
                        ),
                    )
                continue
            candidates.append(
                ShRNACandidate(
                    candidate_id=f"manual-{index}",
                    target_sequence=sequence,
                    start_position=None,
                    intrinsic_score=Decimal("0"),
                    source_rank=index,
                    blast_status=BlastScreenStatus.MANUALLY_ACCEPTED,
                    blast_note="用户手动粘贴并确认",
                ),
            )
        return NewShRNAProjectCommand(
            project_id=self.project_id.text().strip(),
            gene_symbol=self.gene_symbol.text().strip(),
            species=self.species.currentText(),
            cds_sequence=self.cds_sequence.toPlainText(),
            candidates=tuple(candidates),
            target_count=self.target_count.value(),
            clones_per_target=self.clones_per_target.value(),
            received_date=self.received_date.date().toPython(),
            due_date=self.due_date.date().toPython(),
            actor=actor,
            vector_sequence_confirmed=self.vector_confirmation.isChecked(),
            transcript_accession=self.transcript_accession.text().strip() or None,
            primer_vendor_name=_vendor_name_from_combo(self.primer_template_combo),
            sequencing_vendor_name=_vendor_name_from_combo(self.sequencing_template_combo),
            primer_template_id=self.primer_template_combo.currentData(),
            sequencing_template_id=self.sequencing_template_combo.currentData(),
        )


class ImportExpressionProtocolDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("导入表达载体 protocol")
        self.resize(720, 620)
        self._vector: VectorRecord | None = None
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.form = form
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        vector_row = QWidget()
        vector_layout = QHBoxLayout(vector_row)
        vector_layout.setContentsMargins(0, 0, 0, 0)
        self.vector_path = QLineEdit()
        self.vector_path.setReadOnly(True)
        browse = QPushButton("选择图谱")
        browse.clicked.connect(self._browse_vector)
        vector_layout.addWidget(self.vector_path, 1)
        vector_layout.addWidget(browse)

        self.display_name = QLineEdit()
        self.protocol_version_id = QLineEdit()
        self.left_boundary = QSpinBox()
        self.left_boundary.setRange(0, 10_000_000)
        self.right_boundary = QSpinBox()
        self.right_boundary.setRange(1, 10_000_000)
        self.left_homology = QLineEdit()
        self.right_homology = QLineEdit()
        self.kozak = QLineEdit("GCCACC")
        self.stop_rule = QComboBox()
        self.stop_rule.addItem("去掉 stop，与 C 端标签融合", "remove_for_c_terminal_fusion")
        self.stop_rule.addItem("保留 stop，不与 C 端标签融合", "preserve")
        self.fusion_name = QLineEdit("3xFLAG")

        form.addRow("载体图谱", vector_row)
        form.addRow("Protocol 显示名", self.display_name)
        form.addRow("Protocol 版本 ID", self.protocol_version_id)
        form.addRow("左插入边界（0-based）", self.left_boundary)
        form.addRow("右插入边界（0-based）", self.right_boundary)
        form.addRow("F 引物载体同源序列", self.left_homology)
        form.addRow("R 引物载体同源序列", self.right_homology)
        form.addRow("Kozak", self.kozak)
        form.addRow("终止密码子规则", self.stop_rule)
        form.addRow("C 端融合名称", self.fusion_name)
        layout.addLayout(form)

        note = QLabel(
            "当前入口用于导入已经确认过插入边界和同源臂的载体。"
            "保存后，新项目可直接重复使用；载体校验值变化时会拒绝套用。",
        )
        note.setWordWrap(True)
        note.setObjectName("mutedLabel")
        self.note = note
        layout.addWidget(note)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("校验并保存")
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("取消")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _browse_vector(self) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self,
            "选择载体图谱",
            "",
            "载体图谱 (*.dna *.gb *.gbk *.fasta *.fa *.fna)",
        )
        if selected:
            try:
                self.load_vector_path(Path(selected))
            except Exception as error:
                QMessageBox.critical(self, "载体读取失败", str(error))

    def load_vector_path(self, path: Path) -> None:
        records = read_sequence_file(Path(path))
        if len(records) != 1:
            raise ValueError("一个载体文件必须只包含一条序列")
        parsed = records[0]
        provisional = VectorRecord.from_sequence(
            vector_record_id="pending",
            structural_display_name=Path(path).stem,
            sequence=parsed.sequence,
            topology=parsed.topology if parsed.topology != "unknown" else "circular",
            local_aliases=(Path(path).stem,),
        )
        self._vector = VectorRecord.from_sequence(
            vector_record_id=f"local-{provisional.normalized_circular_sha256[:16]}",
            structural_display_name=provisional.structural_display_name,
            sequence=provisional.sequence,
            topology=provisional.topology,
            local_aliases=provisional.local_aliases,
        )
        self.vector_path.setText(str(Path(path)))
        if not self.display_name.text().strip():
            self.display_name.setText(f"{Path(path).stem} 表达 protocol")
        if not self.protocol_version_id.text().strip():
            self.protocol_version_id.setText(
                f"{Path(path).stem}-expression-v1".replace(" ", "-"),
            )

    def profile(self) -> tuple[VectorRecord, ExpressionVectorProtocol]:
        if self._vector is None:
            raise ValueError("请先选择载体图谱")
        display_name = self.display_name.text().strip()
        version_id = self.protocol_version_id.text().strip()
        if not display_name or not version_id:
            raise ValueError("Protocol 显示名和版本 ID 不能为空")
        fusion = self.fusion_name.text().strip() or None
        if self.stop_rule.currentData() == "preserve":
            fusion = None
        protocol = ExpressionVectorProtocol(
            protocol_id=version_id.rsplit("-v", 1)[0],
            protocol_version_id=version_id,
            display_name=display_name,
            status="enabled",
            experimental_validation_status="unverified",
            vector_record_id=self._vector.vector_record_id,
            vector_checksum=self._vector.normalized_circular_sha256,
            workflow_type="expression",
            insertion_mode="confirmed_interval_with_homology_prefixes",
            left_boundary=self.left_boundary.value(),
            right_boundary=self.right_boundary.value(),
            left_primer_homology=self.left_homology.text(),
            right_primer_homology=self.right_homology.text(),
            kozak_sequence=self.kozak.text(),
            stop_codon_rule=self.stop_rule.currentData(),
            c_terminal_fusion_name=fusion,
        )
        validation = validate_expression_protocol(self._vector, protocol)
        if not validation.is_valid:
            raise ValueError("；".join(item.message for item in validation.errors))
        return self._vector, protocol

    def _accept_validated(self) -> None:
        try:
            self.profile()
        except Exception as error:
            QMessageBox.warning(self, "Protocol 尚未通过校验", str(error))
            return
        self.accept()


class NewExpressionProjectDialog(QDialog):
    def __init__(
        self,
        calendar: ChinaBusinessCalendar,
        service: GeneSnapApplicationService,
        parent: QWidget | None = None,
        *,
        ncbi_client: NCBITranscriptClient | None = None,
    ) -> None:
        super().__init__(parent)
        self.calendar = calendar
        self.service = service
        self.ncbi_client = ncbi_client or NCBITranscriptClient()
        self.setWindowTitle("新建表达类项目")
        layout, content_layout = _project_intake_layout(
            self,
            preferred_width=720,
            preferred_height=760,
        )
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self.project_id = QLineEdit()
        self.gene_symbol = QLineEdit()
        self.species = QComboBox()
        self.species.addItems(("human", "mouse", "rat"))
        self.transcript_accession = QLineEdit()
        transcript_row = QWidget()
        transcript_layout = QHBoxLayout(transcript_row)
        transcript_layout.setContentsMargins(0, 0, 0, 0)
        transcript_layout.addWidget(self.transcript_accession, 1)
        self.lookup_transcript_button = QPushButton("查询 CDS")
        self.lookup_transcript_button.clicked.connect(self._lookup_transcript)
        transcript_layout.addWidget(self.lookup_transcript_button)
        self.profile_combo = QComboBox()
        profiles = service.list_expression_profiles()
        if profiles:
            for profile in profiles:
                self.profile_combo.addItem(
                    f"{profile.display_name} · {profile.vector_name}",
                    profile.profile_id,
                )
        else:
            self.profile_combo.addItem("尚未导入表达载体 protocol", None)
        self.cds_sequence = QPlainTextEdit()
        self.cds_sequence.setMinimumHeight(160)
        self.cds_sequence.setPlaceholderText("CDS 或单条 FASTA")
        self.construct_lines = QPlainTextEdit()
        self.construct_lines.setMinimumHeight(120)
        self.construct_lines.setPlaceholderText("一行一个，例如：\nFL\n1-300aa\nΔ301-600\nK436R")
        self.clones_per_construct = QSpinBox()
        self.clones_per_construct.setRange(1, 96)
        self.clones_per_construct.setValue(5)
        self.sequencing_method = QComboBox()
        self.sequencing_method.addItems(("Nanopore", "Sanger"))
        self.primer_template_combo = QComboBox()
        _populate_workbook_template_combo(
            self.primer_template_combo,
            service,
            "primer_order",
        )
        self.sequencing_template_combo = QComboBox()
        _populate_workbook_template_combo(
            self.sequencing_template_combo,
            service,
            "sequencing_order",
        )
        self.received_date = QDateEdit(QDate.currentDate())
        self.due_date = QDateEdit()
        _configure_project_date_edits(self.received_date, self.due_date)
        self.received_date.dateChanged.connect(self._update_due_date)
        self.confirmation_reason = QPlainTextEdit()
        self.confirmation_reason.setMaximumHeight(76)
        self.confirmation_reason.setPlaceholderText(
            "点突变或其他需人工复核的设计，请填写最终序列核对结论",
        )

        form.addRow("项目号", self.project_id)
        form.addRow("基因", self.gene_symbol)
        form.addRow("物种", self.species)
        form.addRow("转录本号（可选）", transcript_row)
        form.addRow("载体 protocol", self.profile_combo)
        form.addRow("CDS / FASTA", self.cds_sequence)
        form.addRow("构建需求", self.construct_lines)
        form.addRow("每个构建送测克隆数", self.clones_per_construct)
        form.addRow("测序方式", self.sequencing_method)
        form.addRow("引物订购模板", self.primer_template_combo)
        form.addRow("送测模板", self.sequencing_template_combo)
        form.addRow("接收日期", self.received_date)
        form.addRow("标准完工日期", self.due_date)
        form.addRow("人工复核说明", self.confirmation_reason)
        content_layout.addLayout(form)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("生成设计并保存")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._update_due_date()

    def _update_due_date(self) -> None:
        _set_default_due_date(
            self.calendar,
            self.received_date,
            self.due_date,
            9,
        )

    def _construct_lines(self) -> tuple[str, ...]:
        return tuple(
            line.strip()
            for line in self.construct_lines.toPlainText().splitlines()
            if line.strip()
        )

    def _lookup_transcript(self) -> None:
        candidate = _lookup_transcript_candidate(
            self,
            self.ncbi_client,
            accession=self.transcript_accession.text(),
            gene_symbol=self.gene_symbol.text(),
            species=self.species.currentText(),
        )
        if candidate is None:
            return
        self.transcript_accession.setText(candidate.accession)
        if candidate.gene_symbol:
            self.gene_symbol.setText(candidate.gene_symbol)
        self.cds_sequence.setPlainText(candidate.cds_sequence)

    def _accept_validated(self) -> None:
        if not self.project_id.text().strip() or not self.gene_symbol.text().strip():
            QMessageBox.warning(self, "缺少项目信息", "项目号和基因不能为空。")
            return
        if self.profile_combo.currentData() is None:
            QMessageBox.warning(self, "缺少载体 protocol", "请先导入表达载体 protocol。")
            return
        if not self.cds_sequence.toPlainText().strip():
            QMessageBox.warning(self, "缺少 CDS", "请粘贴 CDS 或 FASTA。")
            return
        if not self._construct_lines():
            QMessageBox.warning(self, "缺少构建需求", "请至少填写一行构建需求。")
            return
        self.accept()

    def command(self, actor: str) -> NewExpressionProjectCommand:
        profile_id = self.profile_combo.currentData()
        if profile_id is None:
            raise ValueError("尚未选择表达载体 protocol")
        vector, protocol = self.service.load_expression_profile(profile_id)
        return NewExpressionProjectCommand(
            project_id=self.project_id.text().strip(),
            gene_symbol=self.gene_symbol.text().strip(),
            species=self.species.currentText(),
            source_cds=self.cds_sequence.toPlainText(),
            construct_lines=self._construct_lines(),
            received_date=self.received_date.date().toPython(),
            due_date=self.due_date.date().toPython(),
            actor=actor,
            vector=vector,
            protocol=protocol,
            clones_per_construct=self.clones_per_construct.value(),
            transcript_accession=self.transcript_accession.text().strip() or None,
            sequencing_method=self.sequencing_method.currentText(),
            design_confirmation_reason=(
                self.confirmation_reason.toPlainText().strip() or None
            ),
            primer_vendor_name=_vendor_name_from_combo(self.primer_template_combo),
            sequencing_vendor_name=_vendor_name_from_combo(self.sequencing_template_combo),
            primer_template_id=self.primer_template_combo.currentData(),
            sequencing_template_id=self.sequencing_template_combo.currentData(),
        )


class ImportReporterProtocolDialog(ImportExpressionProtocolDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("导入 GL002 reporter 载体 protocol")
        for widget in (self.kozak, self.stop_rule, self.fusion_name):
            label = self.form.labelForField(widget)
            if label is not None:
                label.hide()
            widget.hide()
        self.note.setText(
            "该入口用于保存已确认插入边界和同源臂的 GL002/reporter 载体。"
            "保存后，新项目可直接选择；载体校验值变化时会拒绝套用。",
        )

    def load_vector_path(self, path: Path) -> None:
        super().load_vector_path(path)
        self.display_name.setText(f"{Path(path).stem} reporter protocol")
        self.protocol_version_id.setText(
            f"{Path(path).stem}-reporter-v1".replace(" ", "-"),
        )

    def profile(self) -> tuple[VectorRecord, ReporterVectorProtocol]:
        if self._vector is None:
            raise ValueError("请先选择 reporter 载体图谱")
        display_name = self.display_name.text().strip()
        version_id = self.protocol_version_id.text().strip()
        if not display_name or not version_id:
            raise ValueError("Protocol 显示名和版本 ID 不能为空")
        protocol = ReporterVectorProtocol(
            protocol_id=version_id.rsplit("-v", 1)[0],
            protocol_version_id=version_id,
            display_name=display_name,
            status="enabled",
            experimental_validation_status="unverified",
            vector_record_id=self._vector.vector_record_id,
            vector_checksum=self._vector.normalized_circular_sha256,
            workflow_type="promoter_luciferase_reporter",
            insertion_mode="confirmed_interval_with_homology_prefixes",
            left_boundary=self.left_boundary.value(),
            right_boundary=self.right_boundary.value(),
            left_primer_homology=self.left_homology.text(),
            right_primer_homology=self.right_homology.text(),
            default_sequencing_method="Nanopore",
        )
        validation = validate_reporter_protocol(self._vector, protocol)
        if not validation.is_valid:
            raise ValueError("；".join(item.message for item in validation.errors))
        return self._vector, protocol


class NewReporterProjectDialog(QDialog):
    def __init__(
        self,
        calendar: ChinaBusinessCalendar,
        service: GeneSnapApplicationService,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.calendar = calendar
        self.service = service
        self.setWindowTitle("新建 GL002 promoter-luciferase 项目")
        layout, content_layout = _project_intake_layout(
            self,
            preferred_width=740,
            preferred_height=800,
        )
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        self.project_id = QLineEdit()
        self.gene_symbol = QLineEdit()
        self.species = QComboBox()
        self.species.addItems(("human", "mouse", "rat"))
        self.transcript_accession = QLineEdit()
        self.profile_combo = QComboBox()
        profiles = service.list_reporter_profiles()
        if profiles:
            for profile in profiles:
                self.profile_combo.addItem(
                    f"{profile.display_name} · {profile.vector_name}",
                    profile.profile_id,
                )
        else:
            self.profile_combo.addItem("尚未导入 GL002 reporter protocol", None)
        self.promoter_sequence = QPlainTextEdit()
        self.promoter_sequence.setMinimumHeight(170)
        self.promoter_sequence.setPlaceholderText("按 5'-3' 粘贴 promoter；序列末端最靠近 TSS")
        self.mutation_definitions = QPlainTextEdit()
        self.mutation_definitions.setMinimumHeight(90)
        self.mutation_definitions.setPlaceholderText(
            "一行一个，例如：\nmut1:101-120=ACGT...\nmut2:501-510=TTTT...",
        )
        self.construct_lines = QPlainTextEdit()
        self.construct_lines.setMinimumHeight(110)
        self.construct_lines.setPlaceholderText(
            "一行一个，例如：\nWT\nP1500\nP1000\nP500\nmut1\nmut1+mut2",
        )
        self.clones_per_construct = QSpinBox()
        self.clones_per_construct.setRange(1, 96)
        self.clones_per_construct.setValue(5)
        self.primer_template_combo = QComboBox()
        _populate_workbook_template_combo(
            self.primer_template_combo,
            service,
            "primer_order",
        )
        self.sequencing_template_combo = QComboBox()
        _populate_workbook_template_combo(
            self.sequencing_template_combo,
            service,
            "sequencing_order",
        )
        self.received_date = QDateEdit(QDate.currentDate())
        self.due_date = QDateEdit()
        _configure_project_date_edits(self.received_date, self.due_date)
        self.received_date.dateChanged.connect(self._update_due_date)
        self.confirmation_reason = QPlainTextEdit()
        self.confirmation_reason.setMaximumHeight(72)
        self.confirmation_reason.setPlaceholderText("突变方案请填写最终替换序列复核结论")

        form.addRow("项目号", self.project_id)
        form.addRow("基因", self.gene_symbol)
        form.addRow("物种", self.species)
        form.addRow("转录本号（可选）", self.transcript_accession)
        form.addRow("GL002 protocol", self.profile_combo)
        form.addRow("Promoter 序列", self.promoter_sequence)
        form.addRow("突变定义", self.mutation_definitions)
        form.addRow("构建需求", self.construct_lines)
        form.addRow("每个构建送测克隆数", self.clones_per_construct)
        form.addRow("引物订购模板", self.primer_template_combo)
        form.addRow("送测模板", self.sequencing_template_combo)
        form.addRow("接收日期", self.received_date)
        form.addRow("标准完工日期", self.due_date)
        form.addRow("人工复核说明", self.confirmation_reason)
        content_layout.addLayout(form)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("生成设计并保存")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._update_due_date()

    def _update_due_date(self) -> None:
        _set_default_due_date(
            self.calendar,
            self.received_date,
            self.due_date,
            9,
        )

    @staticmethod
    def _lines(widget: QPlainTextEdit) -> tuple[str, ...]:
        return tuple(
            line.strip()
            for line in widget.toPlainText().splitlines()
            if line.strip()
        )

    def _accept_validated(self) -> None:
        if not self.project_id.text().strip() or not self.gene_symbol.text().strip():
            QMessageBox.warning(self, "缺少项目信息", "项目号和基因不能为空。")
            return
        if self.profile_combo.currentData() is None:
            QMessageBox.warning(self, "缺少 GL002 protocol", "请先导入 GL002 reporter 载体。")
            return
        if not self.promoter_sequence.toPlainText().strip():
            QMessageBox.warning(self, "缺少 promoter", "请粘贴 promoter 序列。")
            return
        if not self._lines(self.construct_lines):
            QMessageBox.warning(self, "缺少构建需求", "请至少填写一行构建需求。")
            return
        self.accept()

    def command(self, actor: str) -> NewReporterProjectCommand:
        profile_id = self.profile_combo.currentData()
        if profile_id is None:
            raise ValueError("尚未选择 GL002 reporter protocol")
        vector, protocol = self.service.load_reporter_profile(profile_id)
        return NewReporterProjectCommand(
            project_id=self.project_id.text().strip(),
            gene_symbol=self.gene_symbol.text().strip(),
            species=self.species.currentText(),
            promoter_sequence=self.promoter_sequence.toPlainText(),
            construct_lines=self._lines(self.construct_lines),
            mutation_definitions=self._lines(self.mutation_definitions),
            received_date=self.received_date.date().toPython(),
            due_date=self.due_date.date().toPython(),
            actor=actor,
            vector=vector,
            protocol=protocol,
            clones_per_construct=self.clones_per_construct.value(),
            transcript_accession=self.transcript_accession.text().strip() or None,
            sequencing_method=protocol.default_sequencing_method,
            design_confirmation_reason=(
                self.confirmation_reason.toPlainText().strip() or None
            ),
            primer_vendor_name=_vendor_name_from_combo(self.primer_template_combo),
            sequencing_vendor_name=_vendor_name_from_combo(self.sequencing_template_combo),
            primer_template_id=self.primer_template_combo.currentData(),
            sequencing_template_id=self.sequencing_template_combo.currentData(),
        )


class DesignConfirmationDialog(QDialog):
    def __init__(
        self,
        prepared: PreparedSYNProject,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("确认 SYN 设计风险")
        self.resize(620, 460)
        layout = QVBoxLayout(self)
        summary = QLabel(
            f"{len(prepared.design.final_sequence)} bp · "
            f"{len(prepared.design.oligos)} 条 oligo · "
            f"{len(prepared.design.module_plan.modules)} 个模块",
        )
        summary.setObjectName("detailTitle")
        layout.addWidget(summary)
        risks = QPlainTextEdit()
        risks.setReadOnly(True)
        risks.setPlainText(
            "\n".join(f"• {item}" for item in prepared.design.design_warnings)
            or "无额外风险项",
        )
        layout.addWidget(risks, 1)
        self.reason = QPlainTextEdit()
        self.reason.setPlaceholderText("人工复核结论或确认原因")
        self.reason.setMaximumHeight(100)
        layout.addWidget(self.reason)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("确认并保存")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_with_reason)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept_with_reason(self) -> None:
        if not self.reason.toPlainText().strip():
            QMessageBox.warning(self, "缺少确认原因", "请填写人工复核结论或原因。")
            return
        self.accept()


_TEMPLATE_FIELD_LABELS = {
    "primer_name": "引物名称",
    "sequence": "引物序列",
    "gene_symbol": "基因名",
    "direction": "方向",
    "length": "长度",
    "purification": "纯化方式",
    "scale": "合成规模",
    "sample_name": "样本名称",
    "clone_no": "克隆号",
    "method": "测序方式",
    "note": "备注",
    "customer_name": "客户姓名",
    "responsible_name": "负责人姓名",
    "organization": "客户单位",
    "phone": "联系电话",
    "email": "邮箱",
    "address": "地址",
    "customer_id": "客户编号",
}


class ImportWorkbookTemplateDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("导入供应商 Excel 模板")
        self.resize(760, 720)
        self._source_path: Path | None = None
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.display_name = QLineEdit()
        self.kind = QComboBox()
        self.kind.addItem("引物订购表", "primer_order")
        self.kind.addItem("测序/送测表", "sequencing_order")
        source_row = QWidget()
        source_layout = QHBoxLayout(source_row)
        source_layout.setContentsMargins(0, 0, 0, 0)
        self.source_path = QLineEdit()
        self.source_path.setReadOnly(True)
        browse = QPushButton("选择 .xlsx")
        browse.clicked.connect(self._browse)
        source_layout.addWidget(self.source_path, 1)
        source_layout.addWidget(browse)
        self.sheet_name = QLineEdit()
        self.header_row = QSpinBox()
        self.header_row.setRange(1, 10000)
        self.data_start_row = QSpinBox()
        self.data_start_row.setRange(1, 10000)
        form.addRow("模板名称", self.display_name)
        form.addRow("模板类型", self.kind)
        form.addRow("模板文件", source_row)
        form.addRow("工作表", self.sheet_name)
        form.addRow("表头行", self.header_row)
        form.addRow("数据起始行", self.data_start_row)
        layout.addLayout(form)
        layout.addWidget(
            QLabel("软件已自动识别数据列；请核对下拉选项，可选“不填写”跳过非必填字段。"),
        )
        self.table_mapping = QTableWidget(0, 2)
        self.table_mapping.setHorizontalHeaderLabels(("软件字段", "写入列"))
        self.table_mapping.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table_mapping.horizontalHeader().setSectionResizeMode(
            0,
            QHeaderView.ResizeMode.ResizeToContents,
        )
        self.table_mapping.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.table_mapping, 1)
        layout.addWidget(
            QLabel("订购信息同样会自动匹配；模板里没有的字段可保持“不填写”。"),
        )
        self.contact_mapping = QTableWidget(0, 2)
        self.contact_mapping.setHorizontalHeaderLabels(("档案字段", "写入位置"))
        self.contact_mapping.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.contact_mapping.horizontalHeader().setSectionResizeMode(
            0,
            QHeaderView.ResizeMode.ResizeToContents,
        )
        self.contact_mapping.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.contact_mapping, 1)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("确认映射并保存")
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("取消")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _browse(self) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self,
            "选择供应商 Excel 模板",
            "",
            "Excel 工作簿 (*.xlsx)",
        )
        if not selected:
            return
        source = Path(selected)
        try:
            inspected = inspect_workbook_template(
                source,
                kind=self.kind.currentData(),
            )
            self._source_path = source
            self._populate_inspection(inspected)
        except Exception as error:
            self._source_path = None
            QMessageBox.warning(self, "模板识别失败", str(error))
            return
        self.source_path.setText(selected)
        if not self.display_name.text().strip():
            self.display_name.setText(source.stem)

    def _populate_inspection(self, inspected: WorkbookTemplateInspection) -> None:
        self.sheet_name.setText(inspected.sheet_name)
        self.header_row.setValue(inspected.header_row)
        self.data_start_row.setValue(inspected.data_start_row)
        if self._source_path is None:
            raise ValueError("请先选择 .xlsx 模板")
        choices = workbook_mapping_choices(
            self._source_path,
            sheet_name=inspected.sheet_name,
            header_row=inspected.header_row,
        )
        table_fields = table_field_names(inspected.kind)
        self.table_mapping.setRowCount(len(table_fields))
        for row, field_name in enumerate(table_fields):
            field_item = QTableWidgetItem(_TEMPLATE_FIELD_LABELS.get(field_name, field_name))
            field_item.setData(Qt.ItemDataRole.UserRole, field_name)
            field_item.setFlags(field_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table_mapping.setItem(row, 0, field_item)
            self.table_mapping.setCellWidget(
                row,
                1,
                self._mapping_combo(
                    choices.table_columns,
                    inspected.table_columns.get(field_name),
                ),
            )
        contact_fields = contact_field_names()
        self.contact_mapping.setRowCount(len(contact_fields))
        for row, field_name in enumerate(contact_fields):
            field_item = QTableWidgetItem(_TEMPLATE_FIELD_LABELS.get(field_name, field_name))
            field_item.setData(Qt.ItemDataRole.UserRole, field_name)
            field_item.setFlags(field_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.contact_mapping.setItem(row, 0, field_item)
            self.contact_mapping.setCellWidget(
                row,
                1,
                self._mapping_combo(
                    choices.contact_cells,
                    inspected.contact_cells.get(field_name),
                ),
            )

    @staticmethod
    def _mapping_combo(
        choices: tuple[tuple[object, str], ...],
        selected_value: object | None,
    ) -> QComboBox:
        combo = QComboBox()
        combo.addItem("不填写", None)
        for value, label in choices:
            combo.addItem(label, value)
        if selected_value not in (None, ""):
            selected_index = combo.findData(selected_value)
            if selected_index < 0:
                combo.addItem(str(selected_value), selected_value)
                selected_index = combo.count() - 1
            combo.setCurrentIndex(selected_index)
        return combo

    @staticmethod
    def _mapping_from_table(table: QTableWidget, *, numeric: bool) -> dict:
        mapping = {}
        for row in range(table.rowCount()):
            field_name = table.item(row, 0).data(Qt.ItemDataRole.UserRole)
            editor = table.cellWidget(row, 1)
            if isinstance(editor, QComboBox):
                value = editor.currentData()
            else:
                value = (table.item(row, 1).text() if table.item(row, 1) else "").strip()
            if value in (None, ""):
                continue
            if numeric:
                column_no = int(value)
                if column_no <= 0:
                    raise ValueError("Excel 列号必须大于 0")
                mapping[field_name] = column_no
            else:
                value = str(value).strip()
                if not re.fullmatch(r"[A-Za-z]{1,3}[1-9][0-9]*", value):
                    raise ValueError(f"无效的 Excel 单元格：{value}")
                mapping[field_name] = value.upper()
        return mapping

    def import_data(self) -> tuple[Path, str, WorkbookTemplateInspection]:
        if self._source_path is None:
            raise ValueError("请先选择 .xlsx 模板")
        display_name = self.display_name.text().strip()
        if not display_name:
            raise ValueError("模板名称不能为空")
        kind = self.kind.currentData()
        table_columns = self._mapping_from_table(self.table_mapping, numeric=True)
        required = {"primer_name", "sequence"} if kind == "primer_order" else {"sample_name"}
        if not required.issubset(table_columns):
            missing = ", ".join(_TEMPLATE_FIELD_LABELS[item] for item in sorted(required - set(table_columns)))
            raise ValueError(f"缺少必需映射：{missing}")
        inspection = WorkbookTemplateInspection(
            kind=kind,
            sheet_name=self.sheet_name.text().strip(),
            header_row=self.header_row.value(),
            data_start_row=self.data_start_row.value(),
            table_columns=table_columns,
            contact_cells=self._mapping_from_table(self.contact_mapping, numeric=False),
        )
        return self._source_path, display_name, inspection

    def _accept_validated(self) -> None:
        try:
            self.import_data()
        except Exception as error:
            QMessageBox.warning(self, "模板映射未完成", str(error))
            return
        self.accept()


class ContactProfileDialog(QDialog):
    def __init__(self, profile: ContactProfile, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("订购信息档案")
        self.resize(560, 430)
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.fields: dict[str, QLineEdit] = {}
        for field_name in contact_field_names():
            editor = QLineEdit(getattr(profile, field_name))
            self.fields[field_name] = editor
            form.addRow(_TEMPLATE_FIELD_LABELS[field_name], editor)
        layout.addLayout(form)
        note = QLabel("所有字段均可留空；只有当前模板中存在并已映射的字段才会写入。")
        note.setObjectName("mutedLabel")
        layout.addWidget(note)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Save,
        )
        buttons.button(QDialogButtonBox.StandardButton.Save).setText("保存档案")
        buttons.button(QDialogButtonBox.StandardButton.Save).setObjectName("primaryButton")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def profile(self) -> ContactProfile:
        return ContactProfile(
            **{field_name: editor.text().strip() for field_name, editor in self.fields.items()},
        )


class DueDateAdjustmentDialog(QDialog):
    def __init__(
        self,
        *,
        received_date: date,
        suggested_due_date: date,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("修正标准完工日期")
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.due_date = QDateEdit(
            QDate(
                suggested_due_date.year,
                suggested_due_date.month,
                suggested_due_date.day,
            ),
        )
        self.due_date.setCalendarPopup(True)
        self.due_date.setDisplayFormat("yyyy-MM-dd")
        self.due_date.setMinimumDate(
            QDate(received_date.year, received_date.month, received_date.day),
        )
        self.note = QPlainTextEdit()
        self.note.setPlaceholderText("例如：修正录入年份；客户调整交付日期")
        self.note.setMaximumHeight(90)
        form.addRow("新完工日期", self.due_date)
        form.addRow("修正原因", self.note)
        layout.addLayout(form)
        hint = QLabel("原始日期会保留在项目记录中，本次修正会写入历史。")
        hint.setObjectName("mutedLabel")
        layout.addWidget(hint)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("确认修正")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self._accept_validated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept_validated(self) -> None:
        if not self.note.toPlainText().strip():
            QMessageBox.warning(self, "缺少修正原因", "请填写标准完工日期的修正原因。")
            return
        self.accept()

    def values(self) -> tuple[date, str]:
        selected = self.due_date.date().toPython()
        return selected, self.note.toPlainText().strip()


class SequencingTrackingDialog(QDialog):
    def __init__(
        self,
        *,
        internal_submission_no: str = "",
        vendor_order_no: str = "",
        correction: bool = False,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("修改送测信息" if correction else "记录送测信息")
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.internal_submission_no = QLineEdit(internal_submission_no)
        self.vendor_order_no = QLineEdit(vendor_order_no)
        self.note = QLineEdit()
        form.addRow("送测编号", self.internal_submission_no)
        form.addRow("订单号", self.vendor_order_no)
        if correction:
            form.addRow("修改说明", self.note)
        layout.addLayout(form)
        hint = QLabel("两项编号均允许暂时留空，收到供应商订单号后可再补录。")
        hint.setObjectName("mutedLabel")
        layout.addWidget(hint)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Ok,
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("保存")
        buttons.button(QDialogButtonBox.StandardButton.Ok).setObjectName("primaryButton")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)


class MainWindow(QMainWindow):
    PROJECT_HEADERS = (
        "项目号",
        "基因/目标",
        "项目大类",
        "状态",
        "送测编号",
        "订单号",
        "接收日期",
        "标准完工日期",
        "剩余工作日",
        "设计摘要",
        "可用克隆",
    )
    DETAIL_TAB_HINTS = (
        "项目汇总由软件自动生成，只读；选择单元格后可按 Ctrl+C 复制。",
        "oligo 是设计结果，只读；选择需要的单元格后可按 Ctrl+C 复制。",
        "实验记录只读：先选择对象，再使用上方可用操作；测序分析会自动写入判定。",
        "历史记录由每次受控操作自动生成，只读且不会覆盖原记录。",
        "双击文件可直接打开；也可以选择单元格后复制文件路径。",
    )

    def __init__(
        self,
        service: GeneSnapApplicationService,
        *,
        actor: str,
        today_provider: Callable[[], date] = date.today,
    ) -> None:
        super().__init__()
        self.service = service
        self.actor = actor
        self.today_provider = today_provider
        self.calendar = ChinaBusinessCalendar.for_2026()
        self.current_project: (
            StoredSYNProject
            | StoredShRNAProject
            | StoredExpressionProject
            | StoredReporterProject
            | None
        ) = None
        self.workflow_service = SYNWorkflowService()
        self.state_service = SYNStateTransitionService()
        self.action_buttons: list[QPushButton] = []
        self.setWindowTitle("GeneSnap Workbench")
        ensure_chinese_font()
        self.resize(1480, 880)
        self.setMinimumSize(1040, 680)
        self.setStyleSheet(APP_STYLESHEET)
        self._build_ui()
        self.refresh_projects()

    def _build_ui(self) -> None:
        toolbar = QToolBar("主工具栏")
        toolbar.setMovable(False)
        toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        self.main_toolbar = toolbar
        self.addToolBar(toolbar)
        title = QLabel("GeneSnap Workbench")
        title.setObjectName("appTitle")
        toolbar.addWidget(title)
        self.visibility_filter = QComboBox()
        self.visibility_filter.addItem("进行中", "active")
        self.visibility_filter.addItem("已完成", "completed")
        self.visibility_filter.addItem("已隐藏", "hidden")
        self.visibility_filter.addItem("全部", "all")
        self.visibility_filter.setToolTip("项目视图")
        self.visibility_filter.currentIndexChanged.connect(self.refresh_projects)
        toolbar.addWidget(self.visibility_filter)
        spacer = QWidget()
        spacer.setSizePolicy(
            spacer.sizePolicy().horizontalPolicy(),
            spacer.sizePolicy().verticalPolicy(),
        )
        toolbar.addSeparator()
        self.new_project_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogNewFolder),
            "新建项目",
            self,
        )
        self.new_project_action.setShortcut(QKeySequence.StandardKey.New)
        self.new_project_action.setToolTip("选择 shRNA、表达类、GL002 或全基因合成")
        self.new_project_action.triggered.connect(self.create_project)
        toolbar.addAction(self.new_project_action)
        self.new_expression_project_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_FileIcon),
            "新建表达项目",
            self,
        )
        self.new_expression_project_action.triggered.connect(
            self.create_expression_project,
        )
        self.import_expression_protocol_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DialogOpenButton),
            "表达载体",
            self,
        )
        self.import_expression_protocol_action.triggered.connect(
            self.import_expression_protocol,
        )
        toolbar.addAction(self.import_expression_protocol_action)
        self.new_reporter_project_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_FileIcon),
            "新建 GL002 项目",
            self,
        )
        self.new_reporter_project_action.triggered.connect(self.create_reporter_project)
        self.import_reporter_protocol_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DialogOpenButton),
            "GL002 载体",
            self,
        )
        self.import_reporter_protocol_action.triggered.connect(
            self.import_reporter_protocol,
        )
        toolbar.addAction(self.import_reporter_protocol_action)
        self.new_syn_project_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_FileIcon),
            "新建 SYN 项目",
            self,
        )
        self.new_syn_project_action.triggered.connect(self.create_syn_project)
        toolbar.addSeparator()
        self.import_workbook_template_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DialogOpenButton),
            "表格模板",
            self,
        )
        self.import_workbook_template_action.triggered.connect(
            self.import_workbook_template,
        )
        toolbar.addAction(self.import_workbook_template_action)
        self.contact_profile_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogInfoView),
            "订购信息",
            self,
        )
        self.contact_profile_action.triggered.connect(self.edit_contact_profile)
        toolbar.addAction(self.contact_profile_action)
        self.save_location_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DirIcon),
            "保存位置",
            self,
        )
        self.save_location_action.setToolTip(str(self.service.projects_root))
        self.save_location_action.triggered.connect(self.change_projects_root)
        toolbar.addAction(self.save_location_action)
        self.refresh_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserReload),
            "刷新",
            self,
        )
        self.refresh_action.setShortcut(QKeySequence.StandardKey.Refresh)
        self.refresh_action.triggered.connect(self.refresh_projects)
        toolbar.addAction(self.refresh_action)
        self.open_folder_action = QAction(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon),
            "打开文件夹",
            self,
        )
        self.open_folder_action.setEnabled(False)
        self.open_folder_action.triggered.connect(self.open_current_folder)
        toolbar.addAction(self.open_folder_action)

        central = QWidget()
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(10, 10, 10, 10)
        root_layout.setSpacing(8)
        self.setCentralWidget(central)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root_layout.addWidget(splitter, 1)
        self.project_table = CopyableTableWidget()
        self.project_table.setColumnCount(len(self.PROJECT_HEADERS))
        self.project_table.setHorizontalHeaderLabels(self.PROJECT_HEADERS)
        self.project_table.setAlternatingRowColors(True)
        self.project_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectItems)
        self.project_table.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection)
        self.project_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.project_table.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOn,
        )
        self.project_table.setHorizontalScrollMode(
            QAbstractItemView.ScrollMode.ScrollPerPixel,
        )
        project_scrollbar = self.project_table.horizontalScrollBar()
        project_scrollbar.setMinimumHeight(20)
        project_scrollbar.setStyleSheet(
            """
            QScrollBar:horizontal {
                height: 20px;
                background: #e7eceb;
                border: 1px solid #c7d0cf;
            }
            QScrollBar::handle:horizontal {
                min-width: 56px;
                margin: 3px;
                border-radius: 4px;
                background: #6f817e;
            }
            QScrollBar::handle:horizontal:hover { background: #526764; }
            QScrollBar::add-page:horizontal,
            QScrollBar::sub-page:horizontal { background: transparent; }
            """,
        )
        self.project_table.verticalHeader().setVisible(False)
        self.project_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents,
        )
        self.project_table.horizontalHeader().setStretchLastSection(True)
        self.project_table.itemSelectionChanged.connect(
            self.on_project_selection_changed,
        )
        splitter.addWidget(self.project_table)

        detail = QWidget()
        detail_layout = QVBoxLayout(detail)
        detail_layout.setContentsMargins(0, 0, 0, 0)
        header = QFrame()
        header.setObjectName("detailHeader")
        header.setFixedHeight(72)
        header_layout = QVBoxLayout(header)
        self.detail_title = QLabel("选择一个项目")
        self.detail_title.setObjectName("detailTitle")
        self.detail_subtitle = QLabel("项目详情和当前操作")
        self.detail_subtitle.setObjectName("mutedLabel")
        header_layout.addWidget(self.detail_title)
        header_layout.addWidget(self.detail_subtitle)
        detail_layout.addWidget(header)
        action_bar = QWidget()
        action_bar.setFixedHeight(42)
        self.action_layout = QHBoxLayout(action_bar)
        self.action_layout.setContentsMargins(0, 4, 0, 4)
        self.action_layout.addStretch(1)
        detail_layout.addWidget(action_bar)
        self.detail_hint = QLabel()
        self.detail_hint.setObjectName("mutedLabel")
        self.detail_hint.setWordWrap(True)
        detail_layout.addWidget(self.detail_hint)
        self.tabs = QTabWidget()
        self.summary_table = CopyableTableWidget(0, 2)
        self.summary_table.setHorizontalHeaderLabels(("字段", "内容"))
        self.summary_table.horizontalHeader().setStretchLastSection(True)
        self.summary_table.verticalHeader().setVisible(False)
        self.summary_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.tabs.addTab(self.summary_table, "项目")
        self.tabs.setTabToolTip(0, self.DETAIL_TAB_HINTS[0])
        self.oligo_table = CopyableTableWidget(0, 6)
        self.oligo_table.setHorizontalHeaderLabels(
            ("名称", "序列（5'-3'）", "方向", "pool", "module", "Tm"),
        )
        self.oligo_table.horizontalHeader().setSectionResizeMode(
            1,
            QHeaderView.ResizeMode.Stretch,
        )
        self.oligo_table.verticalHeader().setVisible(False)
        self.oligo_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.tabs.addTab(self.oligo_table, "oligo")
        self.tabs.setTabToolTip(1, self.DETAIL_TAB_HINTS[1])
        self.experiment_table = CopyableTableWidget(0, 5)
        self.experiment_table.setHorizontalHeaderLabels(
            ("对象", "轮次/步骤", "菌落 PCR", "抽提", "测序确认"),
        )
        self.experiment_table.horizontalHeader().setStretchLastSection(True)
        self.experiment_table.verticalHeader().setVisible(False)
        self.experiment_table.setSelectionBehavior(
            QTableWidget.SelectionBehavior.SelectItems,
        )
        self.experiment_table.setSelectionMode(
            QTableWidget.SelectionMode.ExtendedSelection,
        )
        self.experiment_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.experiment_table.setToolTip(self.DETAIL_TAB_HINTS[2])
        self.experiment_table.cellDoubleClicked.connect(
            self._show_experiment_interaction_hint,
        )
        self.tabs.addTab(self.experiment_table, "实验")
        self.tabs.setTabToolTip(2, self.DETAIL_TAB_HINTS[2])
        self.history_table = CopyableTableWidget(0, 5)
        self.history_table.setHorizontalHeaderLabels(
            ("时间", "事件", "原状态", "新状态", "备注"),
        )
        self.history_table.horizontalHeader().setStretchLastSection(True)
        self.history_table.verticalHeader().setVisible(False)
        self.history_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.tabs.addTab(self.history_table, "历史")
        self.tabs.setTabToolTip(3, self.DETAIL_TAB_HINTS[3])
        self.artifact_table = CopyableTableWidget(0, 4)
        self.artifact_table.setHorizontalHeaderLabels(
            ("类型", "设计版本", "生成时间", "文件"),
        )
        self.artifact_table.horizontalHeader().setStretchLastSection(True)
        self.artifact_table.verticalHeader().setVisible(False)
        self.artifact_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.artifact_table.cellDoubleClicked.connect(self.open_artifact)
        self.tabs.addTab(self.artifact_table, "文件")
        self.tabs.setTabToolTip(4, self.DETAIL_TAB_HINTS[4])
        self.tabs.currentChanged.connect(self._update_detail_hint)
        self._update_detail_hint(self.tabs.currentIndex())
        detail_layout.addWidget(self.tabs, 1)
        splitter.addWidget(detail)
        splitter.setSizes((820, 620))

    def _update_detail_hint(self, index: int) -> None:
        if 0 <= index < len(self.DETAIL_TAB_HINTS):
            self.detail_hint.setText(self.DETAIL_TAB_HINTS[index])

    def _show_experiment_interaction_hint(self, _row: int, _column: int) -> None:
        self.statusBar().showMessage(self.DETAIL_TAB_HINTS[2], 8000)

    def refresh_projects(self) -> None:
        selected_id = self.current_project.project_id if self.current_project else None
        summaries = self.service.list_all_projects()
        visibility_mode = self.visibility_filter.currentData()
        if visibility_mode == "active":
            summaries = tuple(
                item
                for item in summaries
                if not item.is_manually_hidden and item.status != "project_completed"
            )
        elif visibility_mode == "completed":
            summaries = tuple(
                item
                for item in summaries
                if not item.is_manually_hidden and item.status == "project_completed"
            )
        elif visibility_mode == "hidden":
            summaries = tuple(item for item in summaries if item.is_manually_hidden)
        self.project_table.setRowCount(len(summaries))
        today = self.today_provider()
        selected_row = None
        for row, summary in enumerate(summaries):
            remaining = (
                summary.frozen_remaining_workdays
                if summary.status == "abnormal_or_paused"
                and summary.frozen_remaining_workdays is not None
                else self.calendar.remaining_workdays(today, summary.due_date)
            )
            values = (
                summary.project_id,
                summary.target_name,
                summary.project_category,
                self._display_status(
                    summary.status,
                    summary.workflow_type,
                    interruption_type=summary.interruption_type,
                ),
                summary.latest_internal_submission_no,
                summary.latest_vendor_order_no,
                summary.received_date.isoformat(),
                summary.due_date.isoformat(),
                str(remaining),
                summary.design_summary,
                ", ".join(summary.usable_clone_names),
            )
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                if column == 0:
                    item.setData(Qt.ItemDataRole.UserRole, summary.workflow_type)
                if column == 8:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.project_table.setItem(row, column, item)
            self._apply_due_color(row, remaining, summary.status)
            if summary.project_id == selected_id:
                selected_row = row
        if selected_row is not None:
            self.project_table.setCurrentCell(selected_row, 0)
            self.on_project_selection_changed()
        elif summaries:
            self.project_table.setCurrentCell(0, 0)
            self.on_project_selection_changed()
        else:
            self.current_project = None
            self._clear_details()

    @staticmethod
    def _display_status(
        status: str,
        workflow_type: str,
        *,
        interruption_type: str | None = None,
    ) -> str:
        if workflow_type == "de_novo_gene_synthesis":
            return display_status_label(status, workflow_type)
        labels = {
            "recorded": "已录入",
            "design_completed": "设计完成/待订购",
            "primers_ordered": "引物已订购",
            "primers_arrived": "引物已到货",
            "cloning_in_progress": "连接转化中",
            "sequencing_in_progress": "送测中",
            "add_on_in_progress": "加测中",
            "rework_in_progress": "重做中",
            "abnormal_or_paused": "异常" if interruption_type == "abnormal" else "暂停",
            "sequencing_pending_analysis": "测序待分析",
            "analysis_completed": "分析完成",
            "plasmid_prep_in_progress": "复抽/质粒抽提中",
            "plasmid_prep_completed": "质粒抽提完成",
            "project_completed": "项目完成",
        }
        return labels.get(status, status)

    def _apply_due_color(self, row: int, remaining: int, status: str) -> None:
        if status == "project_completed":
            color = QColor("#e5f2ed")
        elif status == "abnormal_or_paused":
            color = QColor("#fff1bd")
        elif remaining < 0:
            color = QColor("#ffd8d4")
        elif remaining <= 3:
            color = QColor("#fff0c7")
        else:
            return
        for column in range(self.project_table.columnCount()):
            self.project_table.item(row, column).setBackground(color)

    def _design_summary(self, stored: StoredSYNProject) -> str:
        route = "模块化" if stored.design.module_plan.route.value == "modular" else "单池"
        return (
            f"SYN：{len(stored.design.final_sequence)} bp / {route} / "
            f"{len(stored.design.oligos)} oligos"
        )

    def _usable_clone_names(self, stored: StoredSYNProject) -> tuple[str, ...]:
        summary = self.workflow_service.get_syn_sequencing_summary(stored.snapshot)
        return summary.usable_clone_names

    def on_project_selection_changed(self) -> None:
        row = self.project_table.currentRow()
        if row < 0 or self.project_table.item(row, 0) is None:
            return
        project_id = self.project_table.item(row, 0).text()
        workflow_type = self.project_table.item(row, 0).data(Qt.ItemDataRole.UserRole)
        self.current_project = self.service.load_any_project(project_id, workflow_type)
        self.open_folder_action.setEnabled(True)
        self._populate_details(self.current_project)

    def _populate_details(
        self,
        stored: (
            StoredSYNProject
            | StoredShRNAProject
            | StoredExpressionProject
            | StoredReporterProject
        ),
    ) -> None:
        if isinstance(stored, StoredReporterProject):
            self._populate_reporter_details(stored)
            return
        if isinstance(stored, StoredExpressionProject):
            self._populate_expression_details(stored)
            return
        if isinstance(stored, StoredShRNAProject):
            self._populate_shrna_details(stored)
            return
        status_label = display_status_label(
            stored.snapshot.status,
            "de_novo_gene_synthesis",
        )
        self.detail_title.setText(f"{stored.project_id} · {stored.target_name}")
        self.detail_subtitle.setText(
            f"{status_label} · {len(stored.design.final_sequence)} bp · "
            f"{stored.design.plasmid_simulation.protocol_version_id}",
        )
        self._populate_summary(stored)
        self._populate_oligos(stored)
        self._populate_experiment(stored)
        self._populate_history(stored)
        self._populate_artifacts(stored)
        self._populate_actions(stored)

    def _populate_reporter_details(self, stored: StoredReporterProject) -> None:
        status_label = self._display_status(
            stored.snapshot.status,
            stored.workflow_type,
            interruption_type=stored.snapshot.interruption_type,
        )
        self.detail_title.setText(f"{stored.project_id} · {stored.gene_symbol}")
        self.detail_subtitle.setText(
            f"{status_label} · {len(stored.design.constructs)} 个构建 · {stored.vector_name}",
        )
        values = (
            ("项目号", stored.project_id),
            ("基因", stored.gene_symbol),
            ("物种", stored.species),
            ("状态", status_label),
            ("接收日期", stored.received_date.isoformat()),
            (
                "标准完工日期",
                (stored.snapshot.effective_due_date or stored.due_date).isoformat(),
            ),
            ("设计版本", stored.design.design_version_id),
            ("GL002 protocol", stored.vector_name),
            ("载体序列校验值", stored.vector_design.vector_checksum),
            ("构建数量", str(len(stored.design.constructs))),
            ("每个构建克隆数", str(stored.snapshot.clones_per_construct)),
            ("项目文件夹", str(stored.project_folder)),
        ) + _sequencing_tracking_rows(stored.snapshot)
        self.summary_table.setRowCount(len(values))
        for row, (label, value) in enumerate(values):
            self.summary_table.setItem(row, 0, QTableWidgetItem(label))
            self.summary_table.setItem(row, 1, QTableWidgetItem(value))

        self.oligo_table.setHorizontalHeaderLabels(
            ("名称", "序列（5'-3'）", "方向", "构建", "退火长度", "载体"),
        )
        plans = {item.construct_id: item for item in stored.vector_design.construct_plans}
        primer_rows = []
        for construct in stored.design.constructs:
            plan = plans[construct.construct_id]
            for primer in (plan.forward_primer, plan.reverse_primer):
                primer_rows.append(
                    (
                        primer.name,
                        primer.sequence,
                        primer.direction,
                        construct.construct_name,
                        str(primer.anneal_length),
                        stored.vector_name,
                    ),
                )
        self.oligo_table.setRowCount(len(primer_rows))
        for row, row_values in enumerate(primer_rows):
            for column, value in enumerate(row_values):
                self.oligo_table.setItem(row, column, QTableWidgetItem(value))

        latest = {}
        for record in stored.snapshot.clone_results:
            latest[record.clone_name] = record
        clone_rows = []
        for construct in stored.design.constructs:
            for clone_name in _project_clone_names(
                stored.snapshot,
                construct.construct_name,
                stored.snapshot.clones_per_construct,
            ):
                result = latest.get(clone_name)
                clone_rows.append(
                    (
                        clone_name,
                        construct.construct_name,
                        result.status.upper() if result is not None else "",
                        "",
                        result.reason if result is not None else "",
                    ),
                )
        self.experiment_table.setRowCount(len(clone_rows))
        for row, row_values in enumerate(clone_rows):
            for column, value in enumerate(row_values):
                self.experiment_table.setItem(row, column, QTableWidgetItem(value))

        history = stored.snapshot.status_history
        self.history_table.setRowCount(len(history))
        for row, event in enumerate(history):
            row_values = (
                event.occurred_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                event.event_type,
                event.from_status or "",
                event.to_status or "",
                event.note or "",
            )
            for column, value in enumerate(row_values):
                self.history_table.setItem(row, column, QTableWidgetItem(value))

        artifacts = self.service.reporter_repository.list_artifacts(stored.project_id)
        self.artifact_table.setRowCount(len(artifacts))
        for row, artifact in enumerate(artifacts):
            row_values = (
                artifact.artifact_type,
                artifact.design_version_id,
                artifact.generated_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                str(artifact.path),
            )
            for column, value in enumerate(row_values):
                item = QTableWidgetItem(value)
                if column == 3:
                    item.setData(Qt.ItemDataRole.UserRole, str(artifact.path))
                self.artifact_table.setItem(row, column, item)
        self._populate_molecular_actions(
            stored,
            analyze_callback=self._analyze_reporter_sequencing,
        )

    def _populate_expression_details(self, stored: StoredExpressionProject) -> None:
        status_label = self._display_status(
            stored.snapshot.status,
            stored.workflow_type,
            interruption_type=stored.snapshot.interruption_type,
        )
        self.detail_title.setText(f"{stored.project_id} · {stored.gene_symbol}")
        self.detail_subtitle.setText(
            f"{status_label} · {len(stored.design.constructs)} 个构建 · "
            f"{stored.vector_name}",
        )
        values = (
            ("项目号", stored.project_id),
            ("基因", stored.gene_symbol),
            ("物种", stored.species),
            ("状态", status_label),
            ("接收日期", stored.received_date.isoformat()),
            (
                "标准完工日期",
                (stored.snapshot.effective_due_date or stored.due_date).isoformat(),
            ),
            ("设计版本", stored.design.design_version_id),
            ("载体 protocol", stored.vector_name),
            ("载体序列校验值", stored.vector_design.vector_checksum),
            ("构建数量", str(len(stored.design.constructs))),
            ("每个构建克隆数", str(stored.snapshot.clones_per_construct)),
            ("项目文件夹", str(stored.project_folder)),
        ) + _sequencing_tracking_rows(stored.snapshot)
        self.summary_table.setRowCount(len(values))
        for row, (label, value) in enumerate(values):
            self.summary_table.setItem(row, 0, QTableWidgetItem(label))
            self.summary_table.setItem(row, 1, QTableWidgetItem(value))
        self.summary_table.resizeColumnToContents(0)

        self.oligo_table.setHorizontalHeaderLabels(
            ("名称", "序列（5'-3'）", "方向", "构建", "退火长度", "载体"),
        )
        plans = {item.construct_id: item for item in stored.vector_design.construct_plans}
        primer_rows = []
        for construct in stored.design.constructs:
            plan = plans[construct.construct_id]
            primer_rows.extend(
                (
                    primer.name,
                    primer.sequence,
                    primer.direction,
                    construct.construct_name,
                    str(primer.anneal_length),
                    stored.vector_name,
                )
                for primer in plan.primers
            )
        self.oligo_table.setRowCount(len(primer_rows))
        for row, values_row in enumerate(primer_rows):
            for column, value in enumerate(values_row):
                self.oligo_table.setItem(row, column, QTableWidgetItem(value))

        latest = {}
        for record in stored.snapshot.clone_results:
            latest[record.clone_name] = record
        clone_rows = []
        for construct in stored.design.constructs:
            for clone_name in _project_clone_names(
                stored.snapshot,
                construct.construct_name,
                stored.snapshot.clones_per_construct,
            ):
                result = latest.get(clone_name)
                result_label = ""
                if result is not None:
                    result_label = result.status.upper()
                    if result.manual_review_status == "usable":
                        result_label += "（人工确认可用）"
                    elif result.manual_review_status == "unusable":
                        result_label += "（人工确认不可用）"
                clone_rows.append(
                    (
                        clone_name,
                        construct.construct_name,
                        result_label,
                        "",
                        result.reason if result is not None else "",
                    ),
                )
        self.experiment_table.setRowCount(len(clone_rows))
        for row, values_row in enumerate(clone_rows):
            for column, value in enumerate(values_row):
                self.experiment_table.setItem(row, column, QTableWidgetItem(value))

        history = stored.snapshot.status_history
        self.history_table.setRowCount(len(history))
        for row, event in enumerate(history):
            history_values = (
                event.occurred_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                event.event_type,
                event.from_status or "",
                event.to_status or "",
                event.note or "",
            )
            for column, value in enumerate(history_values):
                self.history_table.setItem(row, column, QTableWidgetItem(value))

        artifacts = self.service.expression_repository.list_artifacts(stored.project_id)
        self.artifact_table.setRowCount(len(artifacts))
        for row, artifact in enumerate(artifacts):
            artifact_values = (
                artifact.artifact_type,
                artifact.design_version_id,
                artifact.generated_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                str(artifact.path),
            )
            for column, value in enumerate(artifact_values):
                item = QTableWidgetItem(value)
                if column == 3:
                    item.setData(Qt.ItemDataRole.UserRole, str(artifact.path))
                self.artifact_table.setItem(row, column, item)
        self._populate_molecular_actions(
            stored,
            analyze_callback=self._analyze_expression_sequencing,
        )
        unresolved = any(
            record.status == "warning" and record.manual_review_status is None
            for record in latest.values()
        )
        if unresolved and stored.snapshot.status == "analysis_completed":
            self._add_action(
                "确认可用",
                lambda: self._review_expression_clone(True),
            )
            self._add_action(
                "确认不可用",
                lambda: self._review_expression_clone(False),
            )

    def _populate_shrna_details(self, stored: StoredShRNAProject) -> None:
        status_label = self._display_status(
            stored.snapshot.status,
            stored.workflow_type,
            interruption_type=stored.snapshot.interruption_type,
        )
        self.detail_title.setText(f"{stored.project_id} · {stored.gene_symbol}")
        self.detail_subtitle.setText(
            f"{status_label} · {stored.design.target_count} targets · "
            f"{stored.design.vector_protocol_version_id}",
        )
        values = (
            ("项目号", stored.project_id),
            ("基因", stored.gene_symbol),
            ("物种", stored.species),
            ("状态", status_label),
            ("接收日期", stored.received_date.isoformat()),
            (
                "标准完工日期",
                (stored.snapshot.effective_due_date or stored.due_date).isoformat(),
            ),
            ("设计版本", stored.design.design_version_id),
            ("载体序列校验值", stored.design.vector_checksum),
            ("Target 数量", str(stored.design.target_count)),
            ("每个 target 克隆数", str(stored.design.clones_per_target)),
            ("项目文件夹", str(stored.project_folder)),
        ) + _sequencing_tracking_rows(stored.snapshot)
        self.summary_table.setRowCount(len(values))
        for row, (label, value) in enumerate(values):
            self.summary_table.setItem(row, 0, QTableWidgetItem(label))
            self.summary_table.setItem(row, 1, QTableWidgetItem(value))
        self.summary_table.resizeColumnToContents(0)

        self.oligo_table.setHorizontalHeaderLabels(
            ("名称", "序列（5'-3'）", "方向", "Target", "BLAST", "得分"),
        )
        rows = []
        for target in stored.design.targets:
            rows.extend(
                (
                    (
                        target.oligos.forward_name,
                        target.oligos.forward_sequence,
                        "F",
                        f"Target {target.target_no}",
                        target.candidate.blast_status.value,
                        str(target.candidate.intrinsic_score),
                    ),
                    (
                        target.oligos.reverse_name,
                        target.oligos.reverse_sequence,
                        "R",
                        f"Target {target.target_no}",
                        target.candidate.blast_status.value,
                        str(target.candidate.intrinsic_score),
                    ),
                ),
            )
        self.oligo_table.setRowCount(len(rows))
        for row, values_row in enumerate(rows):
            for column, value in enumerate(values_row):
                self.oligo_table.setItem(row, column, QTableWidgetItem(value))

        latest = {}
        for record in stored.snapshot.clone_results:
            latest[record.clone_name] = record
        clone_rows = []
        for target in stored.design.targets:
            prefix = target.clone_names[0].rsplit("-", 1)[0]
            clone_names = set(_project_clone_names(
                stored.snapshot,
                prefix,
                stored.design.clones_per_target,
            ))
            rework_pattern = re.compile(
                rf"^{re.escape(stored.design.gene_symbol)}_"
                rf"{target.target_no}n(?:\d+)?_\d+$",
            )
            clone_names.update(
                sample_name
                for submission in stored.snapshot.sequencing_submissions
                for sample_name in submission.sample_names
                if rework_pattern.fullmatch(sample_name)
            )
            for clone_name in sorted(
                clone_names,
                key=lambda value: (int(re.search(r"(\d+)$", value).group(1)), value),
            ):
                result = latest.get(clone_name)
                clone_rows.append(
                    (
                        clone_name,
                        f"Target {target.target_no}",
                        result.status if result is not None else "",
                        "",
                        result.reason if result is not None else "",
                    ),
                )
        self.experiment_table.setRowCount(len(clone_rows))
        for row, values_row in enumerate(clone_rows):
            for column, value in enumerate(values_row):
                self.experiment_table.setItem(row, column, QTableWidgetItem(value))

        history = stored.snapshot.status_history
        self.history_table.setRowCount(len(history))
        for row, event in enumerate(history):
            history_values = (
                event.occurred_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                event.event_type,
                event.from_status or "",
                event.to_status or "",
                event.note or "",
            )
            for column, value in enumerate(history_values):
                self.history_table.setItem(row, column, QTableWidgetItem(value))

        artifacts = self.service.shrna_repository.list_artifacts(stored.project_id)
        self.artifact_table.setRowCount(len(artifacts))
        for row, artifact in enumerate(artifacts):
            artifact_values = (
                artifact.artifact_type,
                artifact.design_version_id,
                artifact.generated_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                str(artifact.path),
            )
            for column, value in enumerate(artifact_values):
                item = QTableWidgetItem(value)
                if column == 3:
                    item.setData(Qt.ItemDataRole.UserRole, str(artifact.path))
                self.artifact_table.setItem(row, column, item)
        self._populate_molecular_actions(
            stored,
            analyze_callback=self._analyze_shrna_sequencing,
        )
        unresolved = any(
            record.status == "warning" and record.manual_review_status is None
            for record in latest.values()
        )
        if unresolved and stored.snapshot.status == "analysis_completed":
            self._add_action("确认可用", lambda: self._review_shrna_clone(True))
            self._add_action("确认不可用", lambda: self._review_shrna_clone(False))

    def _populate_summary(self, stored: StoredSYNProject) -> None:
        sequencing = self.workflow_service.get_syn_sequencing_summary(
            stored.snapshot,
        )
        values = (
            ("项目号", stored.project_id),
            ("目标名称", stored.target_name),
            ("状态", display_status_label(stored.snapshot.status, "de_novo_gene_synthesis")),
            ("接收日期", stored.received_date.isoformat()),
            ("标准完工日期", stored.due_date.isoformat()),
            ("设计版本", stored.design.design_version_id),
            ("目标序列校验值", stored.design.final_checksum),
            ("路线", stored.design.module_plan.route.value),
            ("模块数", str(len(stored.design.module_plan.modules))),
            ("oligo 数", str(len(stored.design.oligos))),
            ("组装轮次", str(stored.snapshot.syn_assembly_round_no)),
            ("当前子步骤", stored.snapshot.syn_assembly_substep or ""),
            ("测序摘要", sequencing.display_text),
            ("项目文件夹", str(stored.project_folder)),
        )
        self.summary_table.setRowCount(len(values))
        for row, (label, value) in enumerate(values):
            self.summary_table.setItem(row, 0, QTableWidgetItem(label))
            self.summary_table.setItem(row, 1, QTableWidgetItem(value))
        self.summary_table.resizeColumnToContents(0)

    def _populate_oligos(self, stored: StoredSYNProject) -> None:
        self.oligo_table.setHorizontalHeaderLabels(
            ("名称", "序列（5'-3'）", "方向", "pool", "module", "Tm"),
        )
        self.oligo_table.setRowCount(len(stored.design.oligos))
        for row, oligo in enumerate(stored.design.oligos):
            values = (
                oligo.name,
                oligo.sequence,
                oligo.strand,
                oligo.pool_id,
                oligo.module_id,
                str(oligo.tm_metadata.tm_celsius),
            )
            for column, value in enumerate(values):
                self.oligo_table.setItem(row, column, QTableWidgetItem(value))

    def _populate_history(self, stored: StoredSYNProject) -> None:
        history = stored.snapshot.status_history
        self.history_table.setRowCount(len(history))
        for row, event in enumerate(history):
            values = (
                event.occurred_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                event.event_type,
                event.from_status or "",
                event.to_status or "",
                event.note or "",
            )
            for column, value in enumerate(values):
                self.history_table.setItem(row, column, QTableWidgetItem(value))

    def _populate_experiment(self, stored: StoredSYNProject) -> None:
        snapshot = stored.snapshot
        latest_colonies = self._latest_by_clone(snapshot.colonies)
        latest_prep = self._latest_by_clone(snapshot.prep_records)
        latest_confirmations = self._latest_by_clone(
            snapshot.sequencing_confirmations,
        )
        if latest_colonies:
            colonies = sorted(
                latest_colonies.values(),
                key=lambda item: (
                    self._attempt_round(snapshot, item.attempt_id),
                    item.clone_no,
                ),
            )
            self.experiment_table.setRowCount(len(colonies))
            for row, colony in enumerate(colonies):
                round_no = self._attempt_round(snapshot, colony.attempt_id)
                prep = latest_prep.get(colony.clone_id)
                confirmation = latest_confirmations.get(colony.clone_id)
                values = (
                    colony.display_name,
                    f"R{round_no}",
                    colony.result.value,
                    prep.status.value if prep is not None else "",
                    confirmation.result.value if confirmation is not None else "",
                )
                for column, value in enumerate(values):
                    item = QTableWidgetItem(value)
                    item.setData(Qt.ItemDataRole.UserRole, colony.clone_id)
                    self.experiment_table.setItem(row, column, item)
            return
        records = snapshot.step_records
        self.experiment_table.setRowCount(len(records))
        for row, record in enumerate(records):
            values = (
                f"第 {record.step_attempt_no} 次",
                record.substep,
                record.result.value,
                "",
                record.note or "",
            )
            for column, value in enumerate(values):
                self.experiment_table.setItem(row, column, QTableWidgetItem(value))

    def _populate_artifacts(self, stored: StoredSYNProject) -> None:
        artifacts = self.service.repository.list_artifacts(stored.project_id)
        self.artifact_table.setRowCount(len(artifacts))
        for row, artifact in enumerate(artifacts):
            values = (
                artifact.artifact_type,
                artifact.design_version_id,
                artifact.generated_at.astimezone().strftime("%Y-%m-%d %H:%M"),
                str(artifact.path),
            )
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                if column == 3:
                    item.setData(Qt.ItemDataRole.UserRole, str(artifact.path))
                self.artifact_table.setItem(row, column, item)

    def _clear_details(self) -> None:
        self.detail_title.setText("选择一个项目")
        self.detail_subtitle.setText("项目详情和当前操作")
        for table in (
            self.summary_table,
            self.oligo_table,
            self.experiment_table,
            self.history_table,
            self.artifact_table,
        ):
            table.setRowCount(0)
        self._clear_actions()
        self.open_folder_action.setEnabled(False)

    def _clear_actions(self) -> None:
        self.action_buttons.clear()
        while self.action_layout.count() > 1:
            item = self.action_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.hide()
                widget.deleteLater()

    def _populate_actions(self, stored: StoredSYNProject) -> None:
        self._clear_actions()
        label = QLabel(display_status_label(stored.snapshot.status, "de_novo_gene_synthesis"))
        label.setObjectName("mutedLabel")
        self.action_layout.insertWidget(0, label)
        status = stored.snapshot.status
        if status == "design_completed":
            self._add_action("标记 oligo 已订购", self._mark_materials_ordered, primary=True)
        elif status == "materials_ordered":
            self._add_action("标记 oligo 已到货", self._mark_materials_arrived, primary=True)
        elif status == "materials_arrived":
            self._add_action("开始合成组装", self._start_assembly, primary=True)
        elif status == "syn_assembly_in_progress":
            if stored.snapshot.syn_assembly_substep == "colony_pcr":
                if not self._active_round_colonies(stored):
                    self._add_action("生成菌落记录", self._create_colonies, primary=True)
                else:
                    self._add_action("阳性", lambda: self._mark_colonies(SYNColonyPCRResult.POSITIVE))
                    self._add_action("阴性", lambda: self._mark_colonies(SYNColonyPCRResult.NEGATIVE))
                    self._add_action(
                        "不确定",
                        lambda: self._mark_colonies(SYNColonyPCRResult.UNCERTAIN),
                    )
                    self._add_action("选择阳性小提", self._select_clones_for_prep, primary=True)
            else:
                self._add_action("记录本次实验", self._record_assembly_step)
                self._add_action("推进下一步", self._advance_substep, primary=True)
        elif status == "plasmid_prep_in_progress":
            self._add_action("标记小提完成", self._complete_selected_prep, primary=True)
        elif status == "awaiting_sequencing_confirmation":
            self._add_action(
                "确认正确",
                lambda: self._confirm_sequencing(SYNSequencingResult.CORRECT),
                primary=True,
            )
            self._add_action(
                "确认错误",
                lambda: self._confirm_sequencing(SYNSequencingResult.INCORRECT),
            )
            self._add_action("不确定", lambda: self._confirm_sequencing(SYNSequencingResult.UNCERTAIN))
            summary = self.workflow_service.get_syn_sequencing_summary(stored.snapshot)
            if summary.correct_count:
                self._add_action("确认项目完成", self._complete_project, primary=True)
            else:
                self._add_action("追加筛选", self._additional_screening)
                self._add_action("重新合成组装", self._restart_assembly)
        self._add_visibility_action(stored.project_id)

    def _populate_molecular_actions(
        self,
        stored: StoredExpressionProject | StoredShRNAProject | StoredReporterProject,
        *,
        analyze_callback: Callable[[], None],
    ) -> None:
        self._clear_actions()
        status_label = self._display_status(
            stored.snapshot.status,
            stored.workflow_type,
            interruption_type=stored.snapshot.interruption_type,
        )
        label = QLabel(status_label)
        label.setObjectName("mutedLabel")
        self.action_layout.insertWidget(0, label)
        if stored.snapshot.status == "abnormal_or_paused":
            self._add_action("恢复项目", self._resume_molecular_project, primary=True)
            self._add_visibility_action(stored.project_id)
            return
        action_by_status = {
            "design_completed": ("标记引物已订购", "mark_primers_ordered"),
            "primers_ordered": ("标记引物已到货", "mark_primers_arrived"),
            "primers_arrived": ("开始连接转化", "start_cloning"),
            "cloning_in_progress": ("标记为已送测", "mark_sent_for_sequencing"),
            "plasmid_prep_in_progress": ("标记抽提完成", "complete_plasmid_prep"),
            "plasmid_prep_completed": ("确认项目完成", "complete_project"),
        }
        if stored.snapshot.status in action_by_status:
            text, action = action_by_status[stored.snapshot.status]
            self._add_action(
                text,
                lambda checked=False, value=action: self._transition_current_molecular(value),
                primary=True,
            )
        elif stored.snapshot.status == "sequencing_in_progress":
            self._add_action("分析测序", analyze_callback, primary=True)
        elif stored.snapshot.status == "add_on_in_progress":
            self._add_action("分析加测结果", analyze_callback, primary=True)
        elif stored.snapshot.status == "analysis_completed":
            self._add_action("重新分析测序", analyze_callback)
            if self._has_usable_clone_for_each_owner(stored):
                self._add_action("开始质粒抽提", self._start_molecular_prep, primary=True)
            else:
                self._add_action("生成加测送测表", self._generate_addon_sequencing, primary=True)
                self._add_action("重新连接/转化", self._start_molecular_rework)
        elif stored.snapshot.status == "rework_in_progress":
            self._add_action("生成重做送测表", self._generate_rework_submission, primary=True)
        if stored.snapshot.sequencing_submissions:
            self._add_action("修改送测编号", self._edit_latest_sequencing_tracking)
        if stored.snapshot.status != "project_completed":
            self._add_action("修正完工日期", self._adjust_current_due_date)
            self._add_action("暂停/异常", self._mark_molecular_interrupted)
        self._add_visibility_action(stored.project_id)

    def _add_visibility_action(self, project_id: str) -> None:
        hidden = self.service.visibility_store.is_hidden(project_id)
        self._add_action(
            "恢复显示" if hidden else "隐藏项目",
            lambda checked=False, value=not hidden: self._set_current_project_hidden(value),
        )

    def _set_current_project_hidden(self, hidden: bool) -> None:
        if self.current_project is None:
            return
        note, accepted = QInputDialog.getMultiLineText(
            self,
            "隐藏项目" if hidden else "恢复显示",
            "原因或说明",
        )
        if not accepted:
            return
        try:
            self.service.set_project_hidden(
                self.current_project.project_id,
                hidden=hidden,
                reason=note,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "项目显示设置失败", str(error))
            return
        self.current_project = None
        self.visibility_filter.setCurrentIndex(2 if hidden else 0)
        self.refresh_projects()

    def _adjust_current_due_date(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        current_due_date = (
            self.current_project.snapshot.effective_due_date
            or self.current_project.due_date
        )
        suggested_due_date = current_due_date
        if current_due_date < self.current_project.received_date:
            suggested_due_date = self.calendar.add_workdays(
                self.current_project.received_date,
                9,
            )
        dialog = DueDateAdjustmentDialog(
            received_date=self.current_project.received_date,
            suggested_due_date=suggested_due_date,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        new_due_date, note = dialog.values()
        try:
            self.service.adjust_molecular_due_date(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                new_due_date=new_due_date,
                note=note,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "修正完工日期失败", str(error))
            return
        self.refresh_projects()

    @staticmethod
    def _latest_clone_results(stored):
        latest = {}
        for record in stored.snapshot.clone_results:
            latest[record.clone_name] = record
        return latest

    def _has_usable_clone_for_each_owner(
        self,
        stored: StoredExpressionProject | StoredShRNAProject | StoredReporterProject,
    ) -> bool:
        latest = self._latest_clone_results(stored)
        usable = tuple(
            record
            for record in latest.values()
            if record.status == "pass" or record.manually_confirmed_usable
        )
        if isinstance(stored, (StoredExpressionProject, StoredReporterProject)):
            return {item.construct_id for item in usable} == {
                item.construct_id for item in stored.design.constructs
            }
        return {item.target_id for item in usable} == {
            item.target_id for item in stored.design.targets
        }

    def _transition_current_molecular(self, action: str) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        internal_submission_no = ""
        vendor_order_no = ""
        if action == "mark_sent_for_sequencing":
            dialog = SequencingTrackingDialog(parent=self)
            if dialog.exec() != QDialog.DialogCode.Accepted:
                return
            internal_submission_no = dialog.internal_submission_no.text().strip()
            vendor_order_no = dialog.vendor_order_no.text().strip()
        try:
            stored = self.service.transition_molecular_project(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                action=action,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
                internal_submission_no=internal_submission_no,
                vendor_order_no=vendor_order_no,
            )
        except Exception as error:
            QMessageBox.critical(self, "项目状态更新失败", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _mark_molecular_interrupted(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        labels = {"暂停": "pause", "异常": "abnormal"}
        selected, accepted = QInputDialog.getItem(
            self,
            "暂停/异常",
            "请选择类型",
            tuple(labels),
            0,
            False,
        )
        if not accepted:
            return
        note, accepted = QInputDialog.getMultiLineText(
            self,
            f"标记为{selected}",
            "原因或说明",
        )
        if not accepted:
            return
        try:
            stored = self.service.mark_molecular_interrupted(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                interruption_type=labels[selected],
                note=note,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "无法暂停项目", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _resume_molecular_project(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        previous = self.current_project.snapshot.interrupted_previous_status
        statuses = {
            "设计完成/待订购": "design_completed",
            "引物已订购": "primers_ordered",
            "引物已到货": "primers_arrived",
            "连接转化中": "cloning_in_progress",
            "送测中": "sequencing_in_progress",
            "加测中": "add_on_in_progress",
            "分析完成": "analysis_completed",
            "重做中": "rework_in_progress",
            "质粒抽提中": "plasmid_prep_in_progress",
            "质粒抽提完成": "plasmid_prep_completed",
        }
        labels = tuple(statuses)
        default_index = next(
            (index for index, label in enumerate(labels) if statuses[label] == previous),
            0,
        )
        selected, accepted = QInputDialog.getItem(
            self,
            "恢复项目",
            "恢复到状态",
            labels,
            default_index,
            False,
        )
        if not accepted:
            return
        note, accepted = QInputDialog.getMultiLineText(
            self,
            "恢复项目",
            "恢复原因或说明",
        )
        if not accepted:
            return
        try:
            stored = self.service.resume_molecular_project(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
                note=note,
                resume_status=statuses[selected],
            )
        except Exception as error:
            QMessageBox.critical(self, "无法恢复项目", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _edit_latest_sequencing_tracking(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ) or not self.current_project.snapshot.sequencing_submissions:
            return
        latest = self.current_project.snapshot.sequencing_submissions[-1]
        dialog = SequencingTrackingDialog(
            internal_submission_no=latest.internal_submission_no,
            vendor_order_no=latest.vendor_order_no,
            correction=True,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        if not dialog.note.text().strip():
            QMessageBox.warning(self, "缺少修改说明", "补录或修改送测编号时请填写说明。")
            return
        try:
            stored = self.service.update_latest_sequencing_tracking(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                internal_submission_no=dialog.internal_submission_no.text(),
                vendor_order_no=dialog.vendor_order_no.text(),
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
                note=dialog.note.text(),
            )
        except Exception as error:
            QMessageBox.critical(self, "送测信息更新失败", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _generate_addon_sequencing(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        try:
            initial = self.service.preview_addon_sequencing(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
            )
        except Exception as error:
            QMessageBox.warning(self, "暂时不能加测", str(error))
            return
        count, accepted = QInputDialog.getInt(
            self,
            "加测克隆数",
            "每个失败 target/构建新增克隆数",
            initial.clones_per_owner,
            1,
            96,
        )
        if not accepted:
            return
        try:
            preview = self.service.preview_addon_sequencing(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                clones_per_owner=count,
            )
        except Exception as error:
            QMessageBox.warning(self, "无法生成加测预览", str(error))
            return
        summary = (
            f"对象：{', '.join(preview.affected_owner_labels)}\n"
            f"新增：{len(preview.sample_names)} 个克隆\n"
            f"编号：{preview.sample_names[0]} 至 {preview.sample_names[-1]}\n\n"
            "确认后才会生成新送测表并记录本轮加测。"
        )
        if QMessageBox.question(self, "确认加测方案", summary) != QMessageBox.StandardButton.Yes:
            return

        profiles = self.service.list_workbook_templates("sequencing_order")
        options = ["内置标准模板", *(item.display_name for item in profiles)]
        selected_template, accepted = QInputDialog.getItem(
            self,
            "选择送测模板",
            "本轮加测使用：",
            options,
            0,
            False,
        )
        if not accepted:
            return
        template_id = None
        vendor_name = "标准"
        if selected_template != "内置标准模板":
            selected_profile = profiles[options.index(selected_template) - 1]
            template_id = selected_profile.template_id
            vendor_name = selected_profile.display_name
            for suffix in ("测序订购表", "测序表", "送测表"):
                if vendor_name.endswith(suffix):
                    vendor_name = vendor_name[: -len(suffix)].strip()

        tracking = SequencingTrackingDialog(parent=self)
        if tracking.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            stored = self.service.confirm_addon_sequencing(
                preview,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
                sequencing_vendor_name=vendor_name,
                sequencing_template_id=template_id,
                internal_submission_no=tracking.internal_submission_no.text(),
                vendor_order_no=tracking.vendor_order_no.text(),
            )
        except Exception as error:
            QMessageBox.critical(self, "生成加测送测表失败", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _start_molecular_rework(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        note, accepted = QInputDialog.getMultiLineText(
            self,
            "重新连接/转化",
            "请说明本次重做原因和对象",
        )
        if not accepted:
            return
        if not note.strip():
            QMessageBox.warning(self, "缺少重做说明", "重新连接/转化必须填写原因。")
            return
        try:
            stored = self.service.start_molecular_rework(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
                note=note,
            )
        except Exception as error:
            QMessageBox.critical(self, "无法开始重做", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _generate_rework_submission(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        count, accepted = QInputDialog.getInt(
            self,
            "重做送测克隆数",
            "每个重做 target/构建送测克隆数",
            5,
            1,
            96,
        )
        if not accepted:
            return
        try:
            preview = self.service.preview_rework_submission(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                clones_per_owner=count,
            )
        except Exception as error:
            QMessageBox.warning(self, "无法生成重做送测预览", str(error))
            return
        summary = (
            f"对象：{', '.join(preview.affected_owner_labels)}\n"
            f"实验轮次：{self.current_project.snapshot.experiment_attempt_no}\n"
            f"样本：{preview.sample_names[0]} 至 {preview.sample_names[-1]}\n\n"
            "确认后生成送测表，主状态将回到送测中。"
        )
        if QMessageBox.question(self, "确认重做送测方案", summary) != QMessageBox.StandardButton.Yes:
            return
        profiles = self.service.list_workbook_templates("sequencing_order")
        options = ["内置标准模板", *(item.display_name for item in profiles)]
        selected_template, accepted = QInputDialog.getItem(
            self,
            "选择送测模板",
            "本轮重做送测使用：",
            options,
            0,
            False,
        )
        if not accepted:
            return
        template_id = None
        vendor_name = "标准"
        if selected_template != "内置标准模板":
            selected_profile = profiles[options.index(selected_template) - 1]
            template_id = selected_profile.template_id
            vendor_name = selected_profile.display_name
            for suffix in ("测序订购表", "测序表", "送测表"):
                if vendor_name.endswith(suffix):
                    vendor_name = vendor_name[: -len(suffix)].strip()
        tracking = SequencingTrackingDialog(parent=self)
        if tracking.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            stored = self.service.confirm_rework_submission(
                preview,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
                sequencing_vendor_name=vendor_name,
                sequencing_template_id=template_id,
                internal_submission_no=tracking.internal_submission_no.text(),
                vendor_order_no=tracking.vendor_order_no.text(),
            )
        except Exception as error:
            QMessageBox.critical(self, "生成重做送测表失败", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _start_molecular_prep(self) -> None:
        if not isinstance(
            self.current_project,
            (StoredExpressionProject, StoredShRNAProject, StoredReporterProject),
        ):
            return
        latest = self._latest_clone_results(self.current_project)
        usable = tuple(
            record
            for record in latest.values()
            if record.status == "pass" or record.manually_confirmed_usable
        )
        selected = []
        seen_owners = set()
        for record in usable:
            owner = (
                record.construct_id
                if isinstance(
                    self.current_project,
                    (StoredExpressionProject, StoredReporterProject),
                )
                else record.target_id
            )
            if owner not in seen_owners:
                selected.append(record.clone_name)
                seen_owners.add(owner)
        try:
            stored = self.service.transition_molecular_project(
                self.current_project.project_id,
                workflow_type=self.current_project.workflow_type,
                action="start_plasmid_prep",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
                selected_clone_names=tuple(selected),
                note=f"选择抽提克隆：{', '.join(selected)}",
            )
        except Exception as error:
            QMessageBox.critical(self, "开始质粒抽提失败", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _add_action(
        self,
        text: str,
        callback: Callable[[], None],
        *,
        primary: bool = False,
    ) -> None:
        button = QPushButton(text)
        button.setFixedHeight(30)
        button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        if primary:
            button.setObjectName("primaryButton")
        button.clicked.connect(callback)
        self.action_layout.insertWidget(self.action_layout.count() - 1, button)
        self.action_buttons.append(button)

    def _analyze_shrna_sequencing(self) -> None:
        if not isinstance(self.current_project, StoredShRNAProject):
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            outcome = self.service.analyze_shrna_sequencing(
                self.current_project.project_id,
                actor=self.actor,
                analyzed_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "测序分析失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.current_project = outcome.project
        self.refresh_projects()
        unmatched = len(outcome.unmatched_files)
        ambiguous = len(outcome.ambiguous_files)
        QMessageBox.information(
            self,
            "测序分析完成",
            f"分析报告：\n{outcome.analysis_report}\n"
            f"未匹配文件：{unmatched}；歧义文件：{ambiguous}",
        )

    def _analyze_expression_sequencing(self) -> None:
        if not isinstance(self.current_project, StoredExpressionProject):
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            outcome = self.service.analyze_expression_sequencing(
                self.current_project.project_id,
                actor=self.actor,
                analyzed_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "表达类测序分析失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.current_project = outcome.project
        self.refresh_projects()
        QMessageBox.information(
            self,
            "表达类测序分析完成",
            f"分析报告：\n{outcome.analysis_report}\n"
            f"未匹配文件：{len(outcome.unmatched_files)}；"
            f"歧义文件：{len(outcome.ambiguous_files)}",
        )

    def _analyze_reporter_sequencing(self) -> None:
        if not isinstance(self.current_project, StoredReporterProject):
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            outcome = self.service.analyze_reporter_sequencing(
                self.current_project.project_id,
                actor=self.actor,
                analyzed_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "GL002 测序分析失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.current_project = outcome.project
        self.refresh_projects()
        QMessageBox.information(
            self,
            "GL002 测序分析完成",
            f"分析报告：\n{outcome.analysis_report}\n"
            f"未匹配文件：{len(outcome.unmatched_files)}；"
            f"歧义文件：{len(outcome.ambiguous_files)}",
        )

    def _review_expression_clone(self, usable: bool) -> None:
        if not isinstance(self.current_project, StoredExpressionProject):
            return
        row = self.experiment_table.currentRow()
        if row < 0 or self.experiment_table.item(row, 0) is None:
            QMessageBox.warning(self, "请选择克隆", "请先在实验表中选择一个 WARNING 克隆。")
            return
        clone_name = self.experiment_table.item(row, 0).text()
        note, accepted = QInputDialog.getMultiLineText(
            self,
            "人工复核表达类克隆",
            f"{clone_name} · {'确认可用' if usable else '确认不可用'}\n复核说明",
        )
        if not accepted:
            return
        if not note.strip():
            QMessageBox.warning(self, "缺少复核说明", "人工复核必须填写说明。")
            return
        try:
            stored = self.service.confirm_expression_clone_review(
                self.current_project.project_id,
                clone_name=clone_name,
                usable=usable,
                note=note,
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "人工复核失败", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def _review_shrna_clone(self, usable: bool) -> None:
        if not isinstance(self.current_project, StoredShRNAProject):
            return
        row = self.experiment_table.currentRow()
        if row < 0 or self.experiment_table.item(row, 0) is None:
            QMessageBox.warning(self, "请选择克隆", "请先在实验表中选择一个 WARNING 克隆。")
            return
        clone_name = self.experiment_table.item(row, 0).text()
        note, accepted = QInputDialog.getMultiLineText(
            self,
            "人工复核 shRNA 克隆",
            f"{clone_name} · {'确认可用' if usable else '确认不可用'}\n复核说明",
        )
        if not accepted:
            return
        if not note.strip():
            QMessageBox.warning(self, "缺少复核说明", "人工复核必须填写说明。")
            return
        try:
            stored = self.service.confirm_shrna_clone_review(
                self.current_project.project_id,
                clone_name=clone_name,
                usable=usable,
                note=note,
                actor=self.actor,
                reviewed_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "人工复核失败", str(error))
            return
        self.current_project = stored
        self.refresh_projects()

    def action_texts(self) -> tuple[str, ...]:
        return tuple(button.text() for button in self.action_buttons)

    def _persist_snapshot(self, updated) -> None:
        if self.current_project is None:
            return
        expected_revision = self.current_project.snapshot.revision
        self.service.repository.save_snapshot(
            self.current_project.project_id,
            updated,
            expected_revision=expected_revision,
            updated_at=datetime.now().astimezone(),
        )
        project_id = self.current_project.project_id
        self.current_project = self.service.load_project(project_id)
        self.refresh_projects()

    def _run_snapshot_action(self, operation: Callable[[], object]) -> None:
        try:
            updated = operation()
            self._persist_snapshot(updated)
        except Exception as error:
            QMessageBox.critical(self, "操作失败", str(error))

    def _mark_materials_ordered(self) -> None:
        snapshot = self.current_project.snapshot
        self._run_snapshot_action(
            lambda: self.state_service.mark_materials_ordered(
                snapshot,
                event_id=f"ordered-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _mark_materials_arrived(self) -> None:
        snapshot = self.current_project.snapshot
        self._run_snapshot_action(
            lambda: self.state_service.mark_materials_arrived(
                snapshot,
                resuspension_complete=False,
                event_id=f"arrived-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _start_assembly(self) -> None:
        reason, accepted = QInputDialog.getMultiLineText(
            self,
            "材料数据确认",
            "实际 nmol、复溶或正式混池数据未在软件中完整记录。\n继续原因：",
        )
        if not accepted or not reason.strip():
            return
        snapshot = self.current_project.snapshot

        def operation():
            started = self.state_service.start_assembly(
                snapshot,
                MaterialReadiness(
                    is_ready=False,
                    can_start_with_override=True,
                    missing_oligo_ids=(),
                    missing_fields=("resuspension_or_mix_record",),
                    errors=(),
                ),
                confirm_missing=True,
                override_reason=reason,
                event_id=f"start-assembly-{uuid4()}",
                override_id=f"material-override-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            )
            return self.workflow_service.begin_initial_attempt(
                started,
                expected_revision=started.revision,
                attempt_id=f"attempt-{uuid4()}",
                event_id=f"begin-attempt-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            )

        self._run_snapshot_action(operation)

    def _record_assembly_step(self) -> None:
        selected, accepted = QInputDialog.getItem(
            self,
            "记录实验结果",
            "结果",
            ("成功", "失败", "部分成功"),
            editable=False,
        )
        if not accepted:
            return
        mapping = {
            "成功": SYNAssemblyAttemptResult.SUCCESS,
            "失败": SYNAssemblyAttemptResult.FAILED,
            "部分成功": SYNAssemblyAttemptResult.PARTIAL,
        }
        note, note_ok = QInputDialog.getMultiLineText(self, "实验备注", "备注（可空）")
        if not note_ok:
            return
        snapshot = self.current_project.snapshot
        self._run_snapshot_action(
            lambda: self.workflow_service.record_assembly_step(
                snapshot,
                result=mapping[selected],
                note=note or None,
                expected_revision=snapshot.revision,
                record_id=f"step-{uuid4()}",
                event_id=f"step-event-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _advance_substep(self) -> None:
        snapshot = self.current_project.snapshot
        order = (
            "assembly_pcr",
            "amplification_pcr",
            "vector_assembly_transformation",
            "colony_pcr",
        )
        current_index = order.index(snapshot.syn_assembly_substep)
        if current_index + 1 >= len(order):
            return
        target = order[current_index + 1]
        self._run_snapshot_action(
            lambda: self.state_service.advance_assembly_substep(
                snapshot,
                to_substep=target,
                event_id=f"advance-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _create_colonies(self) -> None:
        stored = self.current_project
        default_count = (
            12
            if stored.design.module_plan.route is SYNRoute.MODULAR
            or any(
                risk.severity == "high_risk"
                for risk in stored.design.qc_result.risks
            )
            else 8
        )
        count, accepted = QInputDialog.getInt(
            self,
            "生成菌落记录",
            "菌落数量",
            default_count,
            1,
            96,
        )
        if not accepted:
            return
        snapshot = stored.snapshot
        self._run_snapshot_action(
            lambda: self.workflow_service.create_colonies_for_active_round(
                snapshot,
                target_name=stored.target_name,
                colony_count=count,
                expected_revision=snapshot.revision,
                event_id=f"colonies-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _selected_clone_ids(self) -> tuple[str, ...]:
        clone_ids = []
        for index in self.experiment_table.selectedIndexes():
            item = self.experiment_table.item(index.row(), 0)
            clone_id = item.data(Qt.ItemDataRole.UserRole) if item else None
            if clone_id and clone_id not in clone_ids:
                clone_ids.append(clone_id)
        return tuple(clone_ids)

    def _mark_colonies(self, result: SYNColonyPCRResult) -> None:
        clone_ids = self._selected_clone_ids()
        if not clone_ids:
            QMessageBox.warning(self, "未选择克隆", "请先选择一个或多个菌落。")
            return
        original = self.current_project.snapshot

        def operation():
            updated = original
            for clone_id in clone_ids:
                updated = self.workflow_service.record_colony_pcr(
                    updated,
                    clone_id=clone_id,
                    result=result,
                    observed_note=None,
                    expected_revision=updated.revision,
                    event_id=f"colony-result-{uuid4()}",
                    actor=self.actor,
                    occurred_at=datetime.now().astimezone(),
                )
            return updated

        self._run_snapshot_action(operation)

    def _select_clones_for_prep(self) -> None:
        selected = self._selected_clone_ids() or None
        snapshot = self.current_project.snapshot
        self._run_snapshot_action(
            lambda: self.workflow_service.select_clones_for_prep(
                snapshot,
                clone_ids=selected,
                expected_revision=snapshot.revision,
                event_id=f"select-prep-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _complete_selected_prep(self) -> None:
        clone_ids = self._selected_clone_ids()
        if not clone_ids:
            QMessageBox.warning(self, "未选择克隆", "请选择需要标记小提完成的克隆。")
            return
        original = self.current_project.snapshot

        def operation():
            updated = original
            for clone_id in clone_ids:
                updated = self.workflow_service.record_plasmid_prep(
                    updated,
                    clone_id=clone_id,
                    status=PlasmidPrepStatus.COMPLETED,
                    expected_revision=updated.revision,
                    event_id=f"prep-{uuid4()}",
                    actor=self.actor,
                    occurred_at=datetime.now().astimezone(),
                )
            return updated

        self._run_snapshot_action(operation)

    def _confirm_sequencing(self, result: SYNSequencingResult) -> None:
        clone_ids = self._selected_clone_ids()
        if not clone_ids:
            QMessageBox.warning(self, "未选择克隆", "请选择需要确认的克隆。")
            return
        original = self.current_project.snapshot

        def operation():
            updated = original
            for clone_id in clone_ids:
                updated = self.workflow_service.confirm_sequencing(
                    updated,
                    clone_id=clone_id,
                    result=result,
                    note=None,
                    expected_revision=updated.revision,
                    confirmation_id=f"seq-confirm-{uuid4()}",
                    event_id=f"seq-event-{uuid4()}",
                    actor=self.actor,
                    occurred_at=datetime.now().astimezone(),
                )
            return updated

        self._run_snapshot_action(operation)

    def _complete_project(self) -> None:
        answer = QMessageBox.question(
            self,
            "确认项目完成",
            "已存在人工确认正确的克隆。确认项目完成？",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        snapshot = self.current_project.snapshot
        self._run_snapshot_action(
            lambda: self.workflow_service.complete_project(
                snapshot,
                user_confirmed=True,
                expected_revision=snapshot.revision,
                event_id=f"complete-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _additional_screening(self) -> None:
        count, accepted = QInputDialog.getInt(
            self,
            "追加筛选",
            "新增菌落数量",
            8,
            1,
            96,
        )
        if not accepted:
            return
        stored = self.current_project
        snapshot = stored.snapshot
        try:
            preview = self.workflow_service.preview_additional_screening(
                snapshot,
                target_name=stored.target_name,
                colony_count=count,
                expected_revision=snapshot.revision,
                preview_id=f"preview-{uuid4()}",
            )
        except Exception as error:
            QMessageBox.critical(self, "无法追加筛选", str(error))
            return
        answer = QMessageBox.question(
            self,
            "确认追加筛选",
            f"将新增 {count} 个菌落：\n{preview.display_names[0]} … {preview.display_names[-1]}",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self._run_snapshot_action(
            lambda: self.workflow_service.confirm_additional_screening(
                snapshot,
                preview,
                expected_revision=snapshot.revision,
                event_id=f"additional-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _restart_assembly(self) -> None:
        labels = {
            "Assembly PCR": "assembly_pcr",
            "Amplification PCR": "amplification_pcr",
            "载体组装/转化": "vector_assembly_transformation",
            "菌落 PCR": "colony_pcr",
        }
        selected, accepted = QInputDialog.getItem(
            self,
            "重新合成组装",
            "从哪一步开始",
            tuple(labels),
            editable=False,
        )
        if not accepted:
            return
        snapshot = self.current_project.snapshot
        self._run_snapshot_action(
            lambda: self.workflow_service.restart_assembly(
                snapshot,
                restart_from_substep=labels[selected],
                expected_revision=snapshot.revision,
                attempt_id=f"attempt-{uuid4()}",
                event_id=f"restart-{uuid4()}",
                actor=self.actor,
                occurred_at=datetime.now().astimezone(),
            ),
        )

    def _active_round_colonies(self, stored: StoredSYNProject) -> tuple[object, ...]:
        active_attempt_ids = {
            attempt.attempt_id
            for attempt in stored.snapshot.attempts
            if attempt.syn_assembly_round_no == stored.snapshot.syn_assembly_round_no
        }
        return tuple(
            colony
            for colony in stored.snapshot.colonies
            if colony.attempt_id in active_attempt_ids
        )

    @staticmethod
    def _latest_by_clone(records):
        latest = {}
        for record in records:
            latest[record.clone_id] = record
        return latest

    @staticmethod
    def _attempt_round(snapshot, attempt_id: str) -> int:
        for attempt in snapshot.attempts:
            if attempt.attempt_id == attempt_id:
                return attempt.syn_assembly_round_no
        return 0

    def create_project(self) -> None:
        dialog = NewProjectDialog(self.service.projects_root, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            selected_root = self.service.set_projects_root(Path(dialog.output_root.text()))
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "无法设置保存位置", str(error))
            return
        self.save_location_action.setToolTip(str(selected_root))
        routes = {
            "shrna": self.create_shrna_project,
            "expression": self.create_expression_project,
            "reporter": self.create_reporter_project,
            "syn": self.create_syn_project,
        }
        routes[dialog.workflow_type]()

    def change_projects_root(self) -> None:
        chosen = QFileDialog.getExistingDirectory(
            self,
            "选择 GeneSnap 项目保存位置",
            str(self.service.projects_root),
        )
        if not chosen:
            return
        try:
            selected_root = self.service.set_projects_root(Path(chosen))
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "无法设置保存位置", str(error))
            return
        self.save_location_action.setToolTip(str(selected_root))
        self.statusBar().showMessage(f"项目将保存到：{selected_root}", 10000)

    def create_shrna_project(self) -> None:
        dialog = NewShRNAProjectDialog(
            self.calendar,
            self,
            service=self.service,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            stored = self.service.create_shrna_project(
                dialog.command(self.actor),
                created_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "shRNA 设计失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.current_project = stored
        self.refresh_projects()
        QMessageBox.information(
            self,
            "项目已创建",
            f"项目已保存到：\n{stored.project_folder}",
        )

    def import_workbook_template(self) -> None:
        dialog = ImportWorkbookTemplateDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            source_path, display_name, inspected = dialog.import_data()
            saved = self.service.save_workbook_template(
                source_path,
                display_name=display_name,
                inspected=inspected,
            )
        except Exception as error:
            QMessageBox.critical(self, "保存模板失败", str(error))
            return
        kind_label = "引物订购" if saved.kind == "primer_order" else "测序/送测"
        QMessageBox.information(
            self,
            "模板已保存",
            f"{saved.display_name}\n类型：{kind_label}\n后续新建项目时可直接选择。",
        )

    def edit_contact_profile(self) -> None:
        dialog = ContactProfileDialog(self.service.load_contact_profile(), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            self.service.save_contact_profile(dialog.profile())
        except Exception as error:
            QMessageBox.critical(self, "保存订购信息失败", str(error))
            return
        QMessageBox.information(self, "订购信息已保存", "后续生成模板时会自动填入已映射字段。")

    def import_expression_protocol(self) -> None:
        dialog = ImportExpressionProtocolDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            vector, protocol = dialog.profile()
            saved = self.service.save_expression_profile(vector, protocol)
        except Exception as error:
            QMessageBox.critical(self, "保存表达载体失败", str(error))
            return
        QMessageBox.information(
            self,
            "表达载体已保存",
            f"{saved.display_name}\n校验值：{saved.vector_checksum[:16]}...",
        )

    def create_expression_project(self) -> None:
        if not self.service.list_expression_profiles():
            QMessageBox.warning(
                self,
                "尚无表达载体 protocol",
                "请先点击“导入表达载体”保存至少一个已确认的载体 protocol。",
            )
            return
        dialog = NewExpressionProjectDialog(self.calendar, self.service, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            stored = self.service.create_expression_project(
                dialog.command(self.actor),
                created_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "表达类设计失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.current_project = stored
        self.refresh_projects()
        QMessageBox.information(
            self,
            "项目已创建",
            f"项目已保存到：\n{stored.project_folder}",
        )

    def import_reporter_protocol(self) -> None:
        dialog = ImportReporterProtocolDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            vector, protocol = dialog.profile()
            saved = self.service.save_reporter_profile(vector, protocol)
        except Exception as error:
            QMessageBox.critical(self, "保存 GL002 protocol 失败", str(error))
            return
        QMessageBox.information(
            self,
            "GL002 protocol 已保存",
            f"{saved.display_name}\n校验值：{saved.vector_checksum[:16]}...",
        )

    def create_reporter_project(self) -> None:
        if not self.service.list_reporter_profiles():
            QMessageBox.warning(
                self,
                "尚无 GL002 reporter protocol",
                "请先点击“导入 GL002 载体”保存已确认的载体 protocol。",
            )
            return
        dialog = NewReporterProjectDialog(self.calendar, self.service, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            stored = self.service.create_reporter_project(
                dialog.command(self.actor),
                created_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "GL002 reporter 设计失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.current_project = stored
        self.refresh_projects()
        QMessageBox.information(
            self,
            "项目已创建",
            f"项目已保存到：\n{stored.project_folder}",
        )

    def create_syn_project(self) -> None:
        dialog = NewSYNProjectDialog(self.calendar, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        command = dialog.command(self.actor)
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            prepared = self.service.prepare_syn_project(
                command,
                created_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "设计失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        reason = None
        if prepared.design.requires_confirmation:
            confirm = DesignConfirmationDialog(prepared, self)
            if confirm.exec() != QDialog.DialogCode.Accepted:
                return
            reason = confirm.reason.toPlainText().strip()
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            stored = self.service.save_prepared_syn_project(
                command,
                prepared,
                design_confirmation_reason=reason,
                created_at=datetime.now().astimezone(),
            )
        except Exception as error:
            QMessageBox.critical(self, "保存失败", str(error))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.current_project = stored
        self.refresh_projects()
        QMessageBox.information(
            self,
            "项目已创建",
            f"项目已保存到：\n{stored.project_folder}",
        )

    def open_current_folder(self) -> None:
        if self.current_project is None:
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.current_project.project_folder)))

    def open_artifact(self, row: int, column: int) -> None:
        del column
        item = self.artifact_table.item(row, 3)
        if item is None:
            return
        path = Path(item.data(Qt.ItemDataRole.UserRole))
        if path.exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
        else:
            QMessageBox.warning(self, "文件不存在", str(path))
