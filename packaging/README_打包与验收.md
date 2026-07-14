# Windows 打包与验收

## 构建门槛

`build_windows.ps1` 只有在全量测试通过后才执行 PyInstaller。成品使用单文件、无控制台窗口模式，目标电脑不需要预装 Python。

首次构建先在仓库根目录安装构建依赖：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[build]"
```

脚本优先使用仓库根目录的 `.venv`；没有本地虚拟环境时，才使用 PATH 中的 `python`。也可以通过 `-Python` 显式指定解释器。

## 自动验收

`build_windows.ps1` 会使用独立临时项目目录启动打包后的 GUI EXE，等待进程结束，并检查报告中的 `ok`。任何一步失败都会停止构建。

需要手工复核时，使用 `Start-Process -Wait` 运行；GUI 子系统程序不能依赖普通命令行调用来等待结束：

```powershell
$root = Join-Path $env:TEMP "genesnap-packaged-smoke"
$arguments = @(
  "--smoke-test",
  "--data-dir", ('"{0}"' -f (Join-Path $root "data")),
  "--smoke-report", ('"{0}"' -f (Join-Path $root "report.json"))
)
$process = Start-Process `
  -FilePath "release/GeneSnapWorkbench.exe" `
  -ArgumentList $arguments `
  -WindowStyle Hidden `
  -Wait `
  -PassThru
$process.ExitCode
Get-Content (Join-Path $root "report.json")
```

通过条件：

- `ok = true`。
- 能从 EXE 内部加载公开 pUC57 和 pLKO.1 starter。
- 能分别创建 shRNA、表达类、GL002 reporter 和 SYN 项目。
- 四种工作流均能生成设计、SQLite 记录和至少 5 个输出文件。
- 关闭仓库并重新打开后，四种项目的状态和设计版本保持一致。
- `release/` 同时包含 `LICENSE`、`NOTICE` 和 `THIRD_PARTY_NOTICES.md`，安装版也必须安装这些文件。

## 图形界面验收

1. 双击 `GeneSnapWorkbench.exe` 能启动，不需要 Python 环境。
2. 首屏为项目总览，不是说明页。
3. 工具栏每个按钮都有可见中文名称；点击“新建项目”先选择 shRNA、表达类、GL002 reporter 或 SYN，不直接跳入 shRNA。
4. 首次启动可选择项目保存根目录；新建项目时仍可修改，项目输出不得只落在隐藏的 AppData 目录。
5. shRNA/表达类可通过基因名或转录本号查询 CDS，也可手工粘贴 FASTA。
6. shRNA 默认执行 Broad GPP 候选设计和 NCBI BLAST，界面逐条显示 Broad 得分、位置、BLAST 状态和 oligo 来源。
7. 可导入载体 protocol、供应商 Excel 模板并保存订购信息档案。
8. 关闭程序后重新打开，项目、模板、载体 protocol 和历史记录仍存在。
9. 表格单元格支持复制；双击文件记录能打开对应文件。
10. 测序结果可从项目 `03_sequencing` 文件夹分析，并支持人工复核、加测、重做和抽提完成。
11. 暂停期间剩余工作日冻结，恢复时顺延；完成归档和手动隐藏使用独立视图。
12. shRNA、表达类、GL002 和 SYN 录入窗在短屏幕或 Windows 缩放下可缩小并纵向滚动，底部操作按钮保持可见；项目总览表底部始终有可拖动的横向滚动条。
13. SnapGene `.dna` 能直接导入；供应商模板使用带列字母和表头名称的下拉映射，不要求用户输入数字列号。
14. 日期显示完整年份，完工日期不能早于接收日期；历史错误日期可通过带原因记录的“修正完工日期”操作修正。
15. 项目、oligo、实验、历史和文件页均显示当前可用交互，实验表双击不应被误解为自由编辑。
16. 表达/reporter 载体导入可在“按酶切位点”和“手动粘贴同源臂”之间切换；真实 `.dna` 文件显示格式、长度和拓扑。
17. 供应商引物模板可处理合并单元格，并写入订购日期、Gene ID、NM 号和扩增产物长度。
18. 主表明确区分内部编号、引物送单号、引物订单号、测序送样号和测序订单号；完成项目的编号修改入口被锁定。
19. 质粒抽提必须由用户为每个 Target/构建选择至少一个可用克隆，历史事件和状态显示中文。
20. 主 splitter 可拖动且两侧不会折叠；动作按钮过多时横向滚动；关闭状态的下拉框、数字框和日期框不响应悬停滚轮。

## 当前版本边界

`0.3.4` 是当前可安装整合版本，窗口、任务栏、EXE、安装器和快捷方式使用同一套多尺寸图标。CKO/Flox 仅保留扩展接口；Broad 或 NCBI BLAST 不可用时，shRNA 会明确要求人工确认，不会伪装成已完成脱靶筛查。实验室载体必须由用户导入精确图谱并确认 protocol，软件通过不等于湿实验验证通过。Broad GPP 当前服务条款面向研究用途，商业使用可能需要另行许可。
