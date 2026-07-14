# GeneSnap Workbench

GeneSnap Workbench 是一个 Windows 优先的分子构建工作台。产品目标是把项目管理、序列设计、订购/实验记录、测序判定和可编辑文件输出放在同一个桌面程序中。

## 当前可运行范围

`0.3.4` 已把以下工作流整合到同一个可打包桌面程序，并为窗口、任务栏、EXE、安装器和快捷方式统一加入 GeneSnap Workbench 图标：

- pLKO/shRNA：NCBI CDS 查询、Broad GPP 候选评分、NCBI BLAST 自动脱靶筛选、Broad Full oligo 校验、订购/送测和测序判读；在线服务不可用时才进入明确标记的人工复核备用流程。
- 表达类：全长 CDS、截短体、缺失体、点突变体、多构建、长 CDS 两片段 PCR 和自定义载体 protocol。
- GL002 promoter-luciferase：WT、逐级删除、区域突变和组合突变。
- SYN / de novo 全基因合成：oligo、模块化组装、实验循环和人工测序确认。
- 供应商 Excel 模板自动猜测字段、人工确认映射、保存复用和单一订购信息档案。
- 多轮送测、订单号、加测、重做、WARNING 人工复核、质粒抽提和项目完成。
- 暂停/异常工作日冻结与顺延，以及进行中、已完成、已隐藏三类独立项目视图。
- SQLite 本地持久化、项目文件夹、JSON、XLSX、DOCX 和 GenBank 输出。
- 项目总览、工作日提醒、可复制表格和上下文动作按钮。
- 带中文名称的主工具栏；统一“新建项目”入口先选择四类工作流，不再默认跳到 shRNA。
- 首次启动和每次新建项目均可选择项目保存位置，默认不再把项目输出藏在 C 盘 AppData。
- 四类项目录入窗按屏幕可用空间调整尺寸并提供纵向滚动；项目总览表始终显示清晰的横向滚动条。
- SnapGene `.dna` 图谱按二进制格式读取；表达和 reporter 载体可按酶切位点自动定位同源臂，也可切换为手动粘贴同源臂。
- 供应商模板映射使用带 Excel 字母和表头名称的下拉选项，并支持合并单元格、订购日期、Gene ID、NM 号和扩增产物长度。
- 项目总览分别显示内部编号、引物送单号、引物订单号、测序送样号和测序订单号。
- 质粒抽提前由用户逐个选择真实可用克隆；历史事件和状态统一使用中文显示。
- 主分栏可拖动，动作按钮过多时横向滚动；组合框、数字框和日期框不会因悬停滚轮而误改值。
- 日期统一显示完整年份，禁止保存早于接收日期的完工日期，并提供带原因历史的完工日期修正入口。
- 项目、oligo、实验、历史和文件页明确区分只读记录、复制操作与文件双击打开。

CKO/Flox 在当前版本只保留工作流扩展接口，不包含正式设计实现。内置 pLKO.1 和 pUC57 使用公开参考序列；LV-037、GL002、pcDNA3.1 改造载体等实验室载体通过本地导入并绑定精确序列校验值，不打包内部序列。

## 本地运行

```powershell
$env:PYTHONPATH='src'
python -m genesnap_workbench.app.main
```

运行测试：

```powershell
$env:PYTHONPATH='src'
$env:QT_QPA_PLATFORM='offscreen'
python -m unittest discover -s tests -v
```

## Windows 打包

```powershell
powershell -ExecutionPolicy Bypass -File packaging/build_windows.ps1
```

成品输出到 `release/GeneSnapWorkbench.exe`。构建脚本会先运行全量测试，再调用 PyInstaller，最后生成 SHA-256 文件。

成品自动验收：

```powershell
release/GeneSnapWorkbench.exe `
  --smoke-test `
  --data-dir release/smoke-data `
  --smoke-report release/smoke-report.json
```

## 数据与隐私边界

- 正式代码和公开测试只使用人工构造序列及公开 pUC57 参考序列。
- 不读取、提交或打包 `E:\桌面`、`private_reference/`、公司载体、客户项目、ABI 文件或内部模板。
- 公开 starter 只用于减少录入步骤；真实设计前必须由用户确认实际载体序列一致。
- 软件规则通过不等于湿实验验证通过，`experimental_validation_status` 与软件启用状态分开记录。
- Broad GPP 当前条款面向研究用途，商业使用可能需要另行许可；软件首次在线设计前要求用户阅读并确认条款。

## 目录

```text
src/genesnap_workbench/
  app/                 桌面界面与应用编排
  domain/              不可变领域模型
  project_workflow/    状态机、实验循环和工作日
  resources/           带来源的公开 starter
  sequence_core/       shRNA、表达、reporter、SYN 序列设计内核
  storage/             SQLite 持久化
  template_engine/     JSON/XLSX/DOCX/GenBank 输出
  vector_library/      载体记录、protocol 与序列校验
tests/                 公开单元、端到端和桌面离屏测试
packaging/             Windows 构建与验收脚本
```
