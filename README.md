# GeneSnap Workbench

GeneSnap Workbench 是一个 Windows 优先的分子构建工作台，用于把项目录入、序列设计、订购表输出、实验记录、测序判读和项目进度管理放在同一个桌面程序中。

> 当前版本：`0.3.3 Beta`
>
> 本软件用于辅助设计和流程管理。软件计算通过不代表实验一定成功，所有正式构建方案仍需由具备相应专业能力的人员复核并进行实验验证。

## 下载与安装

仓库发布到 GitHub 后，请从 **Releases** 页面下载 Windows 安装版或便携版。GitHub 页面中的 `Code -> Download ZIP` 只包含源代码，不是可直接使用的安装程序。

- 安装版：运行 `GeneSnapWorkbench_Setup_<版本号>.exe`。
- 便携版：解压后运行 `GeneSnapWorkbench_Portable_<版本号>.exe`。
- 打包版本已经包含 Python 和运行依赖，使用者不需要另外安装 BioPython、PySide6 或 openpyxl。
- Broad GPP、NCBI 转录本查询和 NCBI BLAST 需要能够访问对应在线服务。

## 当前可运行范围

`0.3.3` 已将以下工作流整合到同一个桌面程序，并为窗口、任务栏、EXE、安装器和快捷方式统一加入 GeneSnap Workbench 图标：

- pLKO/shRNA：NCBI CDS 查询、Broad GPP 候选评分、NCBI BLAST 脱靶筛选、oligo 生成、订购/送测和测序判读。
- 表达类：全长 CDS、截短体、缺失体、点突变体、多构建、长 CDS 多片段 PCR 和用户导入载体 protocol。
- promoter-luciferase：WT、逐级删除、区域突变和组合突变。
- SYN / de novo 全基因合成：序列 QC、overlapping oligo、模块化组装、实验循环和人工测序确认。
- 供应商 Excel 模板：自动猜测字段、人工确认映射、保存复用和单一订购信息档案。
- 项目管理：多轮送测、加测、重做、WARNING 人工复核、质粒抽提、暂停/异常工作日顺延和项目完成。
- 本地数据：SQLite 持久化、项目文件夹以及 JSON、XLSX、DOCX 和 GenBank 输出。

CKO/Flox 当前只保留工作流扩展接口，不包含可用的设计页面或正式 protocol。

## 载体与数据边界

- 内置 pLKO.1 和 pUC57 是带来源说明的公开参考序列，真实设计前必须确认其与实际载体一致。
- 实验室改造载体应由用户在本地导入，并绑定精确序列校验值；仓库不提供内部载体序列。
- 仓库及公开测试仅使用公开参考资料或人工构造数据，不包含客户项目、真实测序文件、内部订购模板或个人档案。
- 用户项目数据默认保存在本机，不会因为使用本软件而自动上传到 GitHub。
- Broad GPP 等第三方在线服务受各自条款约束，商业使用前应单独核对当前许可条件。

## 从源码运行

需要 Python 3.11 或更高版本：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e .
$env:PYTHONPATH='src'
.\.venv\Scripts\python.exe -m genesnap_workbench.app.main
```

运行完整测试：

```powershell
$env:PYTHONPATH='src'
$env:QT_QPA_PLATFORM='offscreen'
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

## Windows 打包

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[build]"
powershell -ExecutionPolicy Bypass -File packaging/build_windows.ps1
```

打包脚本会先运行测试，再调用 PyInstaller，并对生成程序执行 smoke test。安装包构建和验收说明见 [packaging/README_打包与验收.md](packaging/README_打包与验收.md)。

## 项目结构

```text
src/genesnap_workbench/
  app/                 桌面界面与应用编排
  domain/              领域模型
  integrations/        Broad、NCBI 等外部服务
  project_workflow/    状态机、实验循环和工作日
  resources/           带来源的公开 starter
  sequence_core/       shRNA、表达、reporter、SYN 序列设计内核
  sequencing/          测序文件匹配与判读
  storage/             SQLite 持久化
  template_engine/     JSON/XLSX/DOCX/GenBank 输出
  vector_library/      载体记录、protocol 与序列校验
tests/                 公开单元、端到端和桌面离屏测试
packaging/             Windows 构建与验收脚本
```

架构说明见 [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)，版本与发布规则见 [docs/发布与版本管理.md](docs/发布与版本管理.md)。

## 贡献与许可

提交修改前请阅读 [CONTRIBUTING.md](CONTRIBUTING.md) 和 [SECURITY.md](SECURITY.md)。本项目代码使用 [Apache License 2.0](LICENSE)；第三方资源和在线服务仍适用其各自的许可与使用条款，详见 [NOTICE](NOTICE) 和 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
