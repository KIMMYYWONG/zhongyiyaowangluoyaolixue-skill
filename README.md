# Network Pharmacology Research Skill

这是一个面向中药网络药理学研究的可审计 Skill。它把网页数据库取数、化合物与靶点标准化、疾病靶点筛选、共病交集、PPI 核心靶点、GO/KEGG、Cytoscape 表和桑基图拆成可独立复用的模块。

核心原则：

- 不把某个方剂、固定药味数量或疾病写死到框架中。
- TCMSP、SYMMAP、HERB2.0 使用不同的来源规则，不混用 OB/DL 门槛。
- 化合物优先按 InChIKey、PubChem CID、CAS 等标识符去重；名称标准化与身份合并分开审计。
- UniProt 映射只接受 Homo sapiens 且 reviewed 的一对一 Gene Symbol；其余进入异常表。
- 疾病数据库先保留原始汇总与来源追踪，再按显式策略生成候选集/推荐集。
- 核心靶点严格按 Betweenness、Closeness、Degree 均大于各自中位数筛选。
- GO 与 KEGG 均使用 clusterProfiler；在线 KEGG 不可用时只能启用显式、带日期的冻结快照模式。
- 完整统计结果、展示子集、QC、参数和会话信息分开输出。

使用时从 [SKILL.md](SKILL.md) 开始。R 脚本位于 `scripts/`，字段模板位于 `assets/`，方法与决策规则位于 `references/`。

真实数据库导出和案例生成物默认不纳入 GitHub。发布前必须确认仓库、可见性、许可证和数据授权。

