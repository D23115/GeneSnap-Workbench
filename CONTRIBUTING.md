# 参与贡献

感谢参与 GeneSnap Workbench。提交代码前，请确保修改可以在不接触任何公司或客户私有数据的情况下复现。

## 开发环境

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e .
```

## 提交前检查

```powershell
$env:PYTHONPATH='src'
$env:QT_QPA_PLATFORM='offscreen'
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

- 新行为应配套测试，修复缺陷应包含能够覆盖原问题的回归测试。
- 不提交真实载体、测序结果、供应商模板、客户项目、个人联系方式或密钥。
- 不把软件计算结果描述为已经通过实验验证。
- 一个 commit 只完成一个清晰目的，推荐使用 `feat:`、`fix:`、`docs:`、`test:`、`refactor:` 或 `chore:` 前缀。

## 问题报告

普通缺陷可以通过 GitHub Issue 报告。安全问题或可能包含敏感数据的问题请按照 [SECURITY.md](SECURITY.md) 私下报告，不要上传原始实验文件。
