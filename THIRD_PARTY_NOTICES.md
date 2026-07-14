# 第三方组件与在线服务说明

GeneSnap Workbench 自有代码使用 Apache License 2.0。下列组件、公开参考数据和在线服务不因进入本项目而改用 Apache-2.0，仍受各自许可证或服务条款约束。

## Python 与打包组件

下表记录 `0.3.2` 本地验收环境中的主要版本。重新发布二进制程序时，应根据实际锁定版本刷新清单，并保留各组件发行包中的完整许可证文件。

| 组件 | 验收版本 | 许可证或许可表达式 | 项目地址 |
| --- | --- | --- | --- |
| BioPython | 1.87 | Biopython License Agreement | https://biopython.org/ |
| docxtpl | 0.20.2 | LGPL-2.1-only | https://github.com/elapouya/python-docx-template |
| openpyxl | 3.1.5 | MIT | https://openpyxl.readthedocs.io/ |
| PySide6 / Qt for Python | 6.11.1 | LGPL-3.0-only OR GPL-2.0-only OR GPL-3.0-only / commercial | https://doc.qt.io/qtforpython-6/ |
| shiboken6 | 6.11.1 | LGPL-3.0-only OR GPL-2.0-only OR GPL-3.0-only | https://doc.qt.io/qtforpython-6/ |
| python-docx | 1.2.0 | MIT | https://python-docx.readthedocs.io/ |
| RapidFuzz | 3.14.5 | MIT | https://github.com/rapidfuzz/RapidFuzz |
| lxml | 6.1.1 | BSD-3-Clause | https://lxml.de/ |
| Jinja2 | 3.1.6 | BSD-3-Clause | https://palletsprojects.com/p/jinja/ |
| MarkupSafe | 3.0.3 | BSD-3-Clause | https://palletsprojects.com/p/markupsafe/ |
| et_xmlfile | 2.0.0 | MIT | https://foss.heptapod.net/openpyxl/et_xmlfile |
| NumPy | 2.5.1 | BSD-3-Clause，并含按其 LICENSE 列出的第三方组件 | https://numpy.org/ |
| PyInstaller | 6.21.0 | GPL-2.0-or-later，带程序分发特别例外 | https://pyinstaller.org/ |

## 公开载体参考

仓库内的 pLKO.1-puro 和 pUC57 GenBank 文件仅保留核苷酸序列、最小 `source` 记录和来源 URL。SnapGene 的图谱、说明和功能注释没有进入仓库。

- pLKO.1-puro 来源页：https://www.snapgene.com/plasmids/viral_expression_and_packaging_vectors/pLKO.1_puro
- pUC57 来源页：https://www.snapgene.com/plasmids/basic_cloning_vectors/pUC57
- SnapGene plasmid resources 条款：https://www.snapgene.com/plasmids

这些参考序列不能替代对实验室实际载体的测序和插入边界确认。

## 在线服务

- Broad GPP：https://portals.broadinstitute.org/gpp/public/
- Broad GPP Terms of Service：https://portals.broadinstitute.org/gppx/portals/public/contactus
- NCBI E-utilities：https://www.ncbi.nlm.nih.gov/books/NBK25501/
- NCBI BLAST：https://blast.ncbi.nlm.nih.gov/

Broad GPP 当前条款将其工具限定为研究用途，并对商业或营利用途提出另行许可要求。软件只提供调用入口和失败后的人工复核路径，使用者必须在实际使用时重新核对最新条款。
