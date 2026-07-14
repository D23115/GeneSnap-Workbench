# AGENTS.md

本文件保存 GeneSnap Workbench 公开仓库的稳定协作规则，供后续 Codex 会话自动读取。

## 工作语言

- 默认用中文与用户沟通。
- README、spec、SOP、计划、ADR 和用户文档必须有中文版本。
- 生物设计自动化与实验验证必须明确区分，不得宣称软件生成方案已经获得实验保证。

## 仓库边界

- 本仓库只允许公开代码、公开来源资源和人工构造测试数据。
- 禁止提交公司载体、客户项目、真实 ABI/AB1、内部 Excel 模板、个人联系方式、API 密钥、本地数据库或带个人信息的截图。
- `.dna`、测序文件、供应商表格和发布二进制默认由 `.gitignore` 排除；确有公开测试需要时，必须先确认来源和脱敏情况，再显式加入。
- 用户导入数据只保存在本机，不得添加任何自动上传功能，除非用户另行明确提出并批准隐私设计。

## 开发流程

- 修改前先阅读 `README.md` 和相关文档、测试。
- 每一轮用户确认的修改完成后，运行与风险相匹配的测试；正式交付前运行完整测试。
- 测试通过并完成隐私检查后，为该轮修改创建一个清晰、单一目的的本地 Git commit。
- 未经用户明确要求，不得执行 `git push`、创建 GitHub Release 或公开仓库。
- 不得把安装包提交到 Git 历史；安装程序、便携版和校验文件应作为 GitHub Release assets 发布。
- 遇到与当前任务无关的本地修改，不得擅自覆盖或回退。

## 产品范围

- 当前实现：shRNA、表达类、promoter-luciferase、SYN/de novo 全基因合成和项目管理。
- CKO/Flox 只保留禁用的扩展接口，不得提供不能实际工作的设计页或占位 protocol。
- 内置公开载体仅作参考；真实设计必须由用户确认实际载体序列和 protocol。

## 验证命令

```powershell
$env:PYTHONPATH='src'
$env:QT_QPA_PLATFORM='offscreen'
python -m unittest discover -s tests -v
```

Windows 发布前还要执行 `packaging/build_windows.ps1` 并运行打包程序 smoke test。
