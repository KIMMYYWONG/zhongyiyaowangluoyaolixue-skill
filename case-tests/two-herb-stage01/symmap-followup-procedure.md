# SYMMAP 后续补录步骤（本轮不执行）

本轮细辛、白芍测试不纳入 SYMMAP 数据。以后补录时严格按以下步骤，不沿用旧表中的 TCMSP MolID 冒充 SYMMAP ingredient ID。

1. 打开 `http://www.symmap.org/search/`，分别搜索细辛和白芍。
2. 核对中药名称、物种或药材实体，记录对应 `SMHB` 编号和详情页 URL；歧义时停止并人工确认。
3. 进入对应 `SMHB` 页面，在 `Related components` 中选择 `ingredient`。
4. 保留完整原始成分表，再按 `OB >= 30` 筛选；记录原始 OB 值、筛选结果和排除原因。
5. 对每个保留成分记录真正的 SYMMAP `INGREDIENT ID`，并进入 ingredient 详情页取得相关 target。
6. 每条记录保存：中药拼音、PubChem CID、SYMMAP ingredient ID、TCMSP MolID（仅在经稳定化学身份确认后填写）、原始成分名、靶点原名、CAS、SMHB ID、来源 URL、访问日期。
7. 以 `InChIKey > 单一 PubChem CID > 无歧义 CAS` 建立跨库化合物身份；多 CID、多 CAS、仅名称记录进入人工审计，不能直接合并。
8. 成分规范名须另经天然产物词典确认；SYMMAP 显示名不能自动标记为 DNP 规范名。
9. 靶点通过 UniProt 限定 `Homo sapiens (9606)` 标准化为 primary Gene Symbol 和 accession；多匹配、非人源、未匹配记录保留在异常表。
10. 经 R 脚本与 TCMSP/HERB 数据合并并去重，同时保留 `source=SYMMAP`、SMHB ID、INGREDIENT ID、原始名称和来源 URL。
