# Full Analysis Walkthrough / 完整分析走讀 (EN / 中文)

---

## 1. Goal 目標

**EN:** Two questions, one AoU carrier cohort (27 cancer-predisposition genes): (1) **Penetrance** — among carriers of a likely-pathogenic variant, what fraction develop the gene's matched cancer, by age? (2) **Concordance** — does ClinVar confirm ACMG-flagged pathogenic variants equally across ancestries?

**中文：** 針對 All of Us 世代中27個癌症易感基因的攜帶者資料，回答兩個問題：(1) **外顯率(Penetrance)** — 帶有「可能致病」變異的人，依年齡層分，有多少比例真的得到該基因對應的癌症？(2) **一致性(Concordance)** — ClinVar 資料庫確認 ACMG 判定為致病的變異，在不同祖源族群中，確認的比例是否一致？

---

## 2. Data 資料

**EN:** Two source tables: `carriers_all_noMUTYH_ancestry.tsv` (one row per person×variant — genotype, ACMG-PLP flag, ClinVar-PLP flag, ancestry) and `All_temp_cohort_simple.csv.gz` (one row per person×condition — age, sex, recorded disease Group; ~306k rows, many people have >1 row because they have multiple recorded conditions).

**中文：** 兩張來源表：`carriers_all_noMUTYH_ancestry.tsv`（每人每個變異一列——genotype、ACMG-PLP標記、ClinVar-PLP標記、祖源）與 `All_temp_cohort_simple.csv.gz`（每人每個疾病紀錄一列——年齡、性別、疾病Group；約30.6萬列，很多人不只一列，因為他們有多筆疾病紀錄）。

---

## 3. Method — Canonical Carriers 方法 — 標準攜帶者名單

**EN:** Every number in this analysis uses ONE fixed carrier list per gene, built the same way every time:
1. Start from every raw variant row for the gene.
2. Collapse to one entry per person per gene.
3. Restrict to **ACMG-PLP** carriers only (not ClinVar-PLP — see rationale below).
4. Join to the cohort table to get age/sex/phenotype; anyone with no matching cohort record is dropped **and that drop is counted, never silent** (roughly 30% per gene).
5. Restrict to the correct sex if the gene's cancer type is sex-specific (e.g. BRCA1 → female only).

This fixed chain matters because earlier, un-canonicalized versions of this pipeline produced three different BRCA1 carrier counts (565 / 611 / 917) in three different tables — each one quietly building its own denominator. Canonicalizing means one number, reused everywhere.

**中文：** 這個分析裡的每一個數字，都用「同一份、每個基因只有一份」的攜帶者名單，用固定流程建立：
1. 從該基因的每一筆原始變異紀錄開始。
2. 合併成「每人每基因」一筆。
3. 只保留 **ACMG-PLP**（不是ClinVar-PLP——原因見下）攜帶者。
4. 跟世代表合併取得年齡/性別/表型；沒有對應世代紀錄的人會被排除，**而且這個排除數字會被明確記錄，不會悄悄消失**（每個基因大約流失30%）。
5. 如果該基因對應的癌症是性別相關，就限制在該性別（例如BRCA1只算女性）。

這個固定流程很重要，因為早期還沒標準化之前，同一個BRCA1在三張不同表格裡出現過三個不同的攜帶者人數（565 / 611 / 917）——每張表都各自算出稍微不同的分母。標準化之後就只有一個數字，到處通用。

---

## 4. Method — Why ACMG, not ClinVar 方法 — 為何選ACMG而非ClinVar

**EN:** ACMG-PLP is a rule-based classification applied the same way regardless of ancestry. ClinVar-PLP depends on manual submissions historically concentrated in European-ancestry variants. Using ClinVar-PLP as the carrier definition would import that ancestry gap directly into the penetrance denominator — which is exactly what Finding 1 (below) measures. So ACMG is used for all headline numbers; a **separate comparison report** (`report_acmg_vs_clinvar.html`) reruns the identical pipeline under ClinVar-PLP to show what that choice would have cost.

**中文：** ACMG-PLP 是規則式分類，不論族群都用同一套規則判定。ClinVar-PLP 依賴人工提交的紀錄，historically集中在歐洲祖源的變異上。如果用ClinVar-PLP當作攜帶者定義，就會把這個祖源落差直接帶進外顯率的分母裡——而這正是下面Finding 1要測量的東西。所以所有headline數字都用ACMG；另外有一份**獨立的比較報告**（`report_acmg_vs_clinvar.html`），把整套流程原封不動地用ClinVar-PLP重跑一次，用來呈現「如果當初選ClinVar會付出什麼代價」。

---

## 5. Method — Penetrance & Background Formulas 方法 — 外顯率與背景比率公式

**EN:**
```
Penetrance(gene, cancer, age band) = N(canonical carriers with that matched cancer) / N(all canonical carriers in that age band)

Background rate(disease, age band) = N(whole cohort with that disease) / N(whole cohort in that age band)
```
Each cancer type is its own line (BRCA1:Breast and BRCA1:Ovary are never merged). A carrier is only "affected" if their recorded disease matches the gene's own cancer-type mapping — a BRCA1 carrier with lung cancer is not a BRCA1 event. Family-history-only records are excluded from the numerator by an assertion in the code. Both proportions get a 95% **Wilson score interval** (chosen over normal approximation because many cells have small N or proportions near 0/1). A cell is marked ★ "significantly elevated" only when the carrier rate's CI lower bound exceeds the background rate's CI upper bound.

**中文：**
```
外顯率(基因, 癌症, 年齡層) = 該年齡層內、有記錄該癌症的標準攜帶者人數 ÷ 該年齡層內全部標準攜帶者人數

背景比率(疾病, 年齡層) = 整個世代中有該疾病的人數 ÷ 整個世代在該年齡層的總人數
```
每種癌症分開算一條線（BRCA1的乳癌跟卵巢癌絕不合併）。只有記錄的疾病跟該基因對應的癌症類型吻合，才算「患病」——BRCA1帶原者得肺癌不算BRCA1事件。只有家族史紀錄（本人未患病）會被程式明確排除在分子之外。兩個比例都用95% **Wilson信賴區間**（比常態近似更適合小樣本或比例接近0/1的情況）。只有當攜帶者比率的信賴區間下界超過背景比率的信賴區間上界時，才標記★「顯著高於背景」。

---

## 6. Method — Age Bins, Artifact Exclusion, Thresholds 方法 — 年齡分組、偽基因排除、門檻

**EN:** Age bins: `<40, 40-49, 50-59, 60-69, 70+` (current age, not age-at-diagnosis; left-closed/right-open, e.g. 50 falls in "50-59"). PMS2 is always excluded (known PMS2/PMS2CL pseudogene mismapping artifact); other genes are flagged the same way if one non-EUR ancestry shows >2× the EUR carrier count. A gene needs **≥20 canonical carriers**, no artifact flag, and a defined cancer mapping to appear in main results — 23 of 27 genes qualify. The same ≥20 bar (applied to both classifiers) gates the ACMG-vs-ClinVar comparison figures — 15 of 27 genes qualify there.

**中文：** 年齡分組：`<40, 40-49, 50-59, 60-69, 70+`（用當前年齡，不是診斷時年齡；左閉右開，例如剛滿50歲會落在「50-59」這組）。PMS2 一律排除（已知的PMS2/PMS2CL偽基因誤定位問題）；其他基因如果有一個非歐裔族群的攜帶者數超過歐裔的兩倍，也會用同樣邏輯被標記。一個基因要出現在主要結果裡，需要 **≥20個標準攜帶者**、沒有被標記為偽基因異常、且有明確的基因-癌症對應——27個基因中有23個符合。同樣的≥20門檻（兩邊classifier都要達標）用來篩選ACMG vs ClinVar比較報告要顯示哪些圖——27個基因中有15個符合。

---

## 7. Main Result 1 — Concordance Gap 主要結果一 — 一致性落差

**EN:** Among ACMG-PLP carriers, the fraction with **no** matching ClinVar-PLP confirmation, by ancestry:

| Ancestry 祖源 | Discordance 不一致率 |
|---|---|
| EUR 歐裔 | 42.5% |
| AFR 非裔 | 51.9% |
| AMR 美洲裔 | 53.2% |
| SAS 南亞裔 | 67.8% |
| EAS 東亞裔 | 90.9% |

Every non-EUR confidence interval sits clearly above EUR's — not noise. **Caveat**: the EAS 90.9% figure is 83.5% attributable to a single recurrent MSH6 variant (chr2:47806842) — read as that one variant's classification status, not broad EAS under-curation. Other ancestries' gaps are spread across many variants — the more typical under-curation pattern.

**中文：** 在ACMG-PLP攜帶者中，「沒有」對應ClinVar-PLP確認的比例，依祖源分：EUR 42.5%、AFR 51.9%、AMR 53.2%、SAS 67.8%、EAS 90.9%。每個非歐裔族群的信賴區間都明顯高於歐裔——不是雜訊。**但要注意**：EAS這個90.9%的數字，有83.5%來自單一一個反覆出現的MSH6變異（chr2:47806842）——應該解讀成「這一個變異的分類狀態」，而不是東亞裔整體被廣泛忽略。其他族群的落差則是分散在很多個不同變異上——這才是比較典型的「資料庫收錄不足」模式。

---

## 8. Main Result 2 — Penetrance 主要結果二 — 外顯率

**EN:** BRCA1: N=320 canonical carriers, **44.7%** develop breast/ovarian cancer, rising with age, above background at every band. BRCA2: N=705, **29.4%**, same shape, consistently lower than BRCA1 (matches known biology). All 23 interpretable genes have their own age-stratified chart. Excluded: BMPR1A (n=13), NF2 (n=19), WT1 (n=16) — below the 20-carrier bar; PMS2 (n=263) — artifact. These numbers run lower than clinic-ascertained estimates (which select for high-risk families already flagged for testing) — expected for an unselected volunteer cohort, not a red flag.

**中文：** BRCA1：320個標準攜帶者，**44.7%**得到乳癌或卵巢癌，隨年齡上升，每個年齡層都高於背景值。BRCA2：705個攜帶者，**29.4%**，型態相同，一直低於BRCA1（符合已知生物學）。23個可解讀的基因都各自有年齡分層的圖表。被排除的：BMPR1A（n=13）、NF2（n=19）、WT1（n=16）——低於20人門檻；PMS2（n=263）——偽基因異常。這些數字比臨床篩選出來的估計值（通常來自已被懷疑高風險而做檢測的家族）來得低——這是預期中的，因為AoU是未經篩選的志願世代，不是警訊。

---

## 9. Main Result 3 — ACMG vs ClinVar Comparison 主要結果三 — ACMG對比ClinVar

**EN:** Rerunning the identical pipeline with ClinVar-PLP instead of ACMG-PLP: total carrier pool drops from **3,790 to 2,066 — a 45.5% loss**, concentrated unevenly. Worst hit: TSC2 keeps only 2.3% (306→7), SMAD4 2.2% (134→3), MAX 5.1% (39→2). Least hit: SDHC 84.5%, SDHB 82.2%, BRCA2 80.0%. **BRCA1 is an exception** — ClinVar actually gives more carriers (381 vs. 320, penetrance 46.7% vs. 44.7%) — this doesn't contradict the ancestry finding, it just means ClinVar happens to classify more BRCA1 variants overall; the ancestry pattern is about *which* carriers get dropped, not every gene's raw total. This turns Finding 1's abstract percentage into a concrete carrier-count cost.

**中文：** 把整套流程原封不動改用ClinVar-PLP重跑一次：攜帶者總數從**3,790降到2,066——流失45.5%**，而且分布很不平均。流失最嚴重：TSC2只剩2.3%（306→7）、SMAD4剩2.2%（134→3）、MAX剩5.1%（39→2）。流失最少：SDHC 84.5%、SDHB 82.2%、BRCA2 80.0%。**BRCA1是例外**——ClinVar反而給出更多攜帶者（381 vs 320，外顯率46.7% vs 44.7%）——這不代表祖源落差的結論有問題，只是ClinVar剛好把更多BRCA1變異歸類為致病；祖源那個發現講的是「哪些攜帶者會被剔除」，不是每個基因的原始總數都會往同一個方向變化。這份比較報告把Finding 1那個抽象的百分比，轉換成具體、可以指出來的攜帶者流失數字。

---

## 10. Limitations 限制

**EN:** Crude, cross-sectional proportions — not survival/cumulative-incidence estimates. AoU is a healthcare-utilizing volunteer cohort, not a random population sample. ~30% of carriers per gene dropped for no matching cohort record (reported, not hidden). Pseudogene-artifact rule is a heuristic, not exhaustive. Single-cohort, preliminary, not externally validated.

**中文：** 這是粗略的橫斷面比例，不是存活分析或累積發生率估計。AoU是有使用醫療資源的志願世代，不是隨機人口抽樣。每個基因大約30%的攜帶者因為沒有對應世代紀錄而被排除（有明確記錄，不是隱藏）。偽基因異常判定規則是經驗法則，不是窮盡式的品管。單一世代、屬於初步結果，尚未經外部資料驗證。
