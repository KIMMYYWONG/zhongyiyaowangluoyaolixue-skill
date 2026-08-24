# 细辛与白芍 Stage 01 TCMSP-only 测试报告

## 本轮范围

- 仅测试细辛与白芍。
- 本轮不抓取 SYMMAP；后续操作已保存在 `symmap-followup-procedure.md`。
- 暂不纳入 HERB2.0 补充关系。
- TCMSP 筛选条件：OB >= 30%，DL >= 0.18。

## 暂定结果

| 药物 | 通过 TCMSP 阈值的来源成分 | 暂定成分-靶点关系 | UniProt 标准靶点 |
|---|---:|---:|---:|
| 白芍 | 8 | 126 | 94 |
| 细辛 | 8 | 173 | 111 |

合计 16 条药物-来源成分记录，对应 15 个唯一化学实体；kaempferol 为两味药共有成分。暂定关系共 299 条，合并去重后为 128 个标准基因，其中两味药共有 77 个。

## 正式冻结状态

当前只能作为 **TCMSP-only 暂定结果**，不能冻结为 Stage 01 正式结果。16 个候选成分均已有稳定化合物主键并通过 TCMSP 阈值，已接受的靶点均使用 reviewed human UniProt 映射；唯一未完成的硬门槛是天然产物词典（DNP）规范名核验。因此正式输出表保持为空，避免把 TCMSP/PubChem 暂用名冒充 DNP 规范名。

## 本轮发现并修正的问题

原规则会把 `DNP review pending` 中出现的 `DNP` 误判为已核验。本轮已增加否定条件：命名依据不得包含 pending、not confirmed、unconfirmed、待核验、待确认、未核验或未确认。

## 三遍验证

1. 真实案例规则检查：16 个成分的身份、阈值和 UniProt 关系均通过；DNP 门槛按预期失败。
2. 标准阳性夹具回归：3 个成分、3 条成分-靶点关系、3 个唯一靶点全部通过，说明修正规则没有破坏合法输入。
3. 完整重复运行：暂定关系表、汇总表及正式 QC 表的 SHA-256 均与首次运行一致。

## 主要输出

- `tcmsp-only-inputs/provisional_tcmsp_compound_target.csv`：可查看的 TCMSP-only 暂定结果。
- `tcmsp-only-inputs/provisional_tcmsp_summary.csv`：两味药汇总。
- `tcmsp-only-stage01/07_stage01_qc.csv`：正式门槛检查。
- `tcmsp_compound_evidence.csv`：16 个来源成分的 OB、DL、稳定标识及来源链接。
- `symmap-followup-procedure.md`：后续 SYMMAP 补充步骤。

