# DNAFlow v0.0.1 - 高效模块化全基因组分析流程

**DNAFlow** 是一款基于 Snakemake 开发的、面向工业级应用的 WGS/WES 自动化分析流程。它遵循 GATK Best Practices 和群体遗传学标准，旨在提供高度灵活、可扩展且透明的基因组变异检测方案。

## 🚀 今日更新摘要 (2026-03-05)

### 1. 变异检测模块化重构 (Germline Variant Calling)
- **多引擎支持**：支持通过 `config.yaml` 动态切换 Bcftools、GATK Single-Sample (GVCF 模式) 和 GATK Joint-Calling。
- **动态过滤算法**：为 Bcftools 引入了基于样本量的动态 DP 阈值（`num_samples * 5`）及 MAF/Missingness 自动切换逻辑。

### 2. Somatic（体细胞）分析链条全线贯通
- **Mutect2 配对分析**：实现了从 `Tumor-Normal` 配对调用到多重过滤（Orientation Bias, Contamination, Segmentation）的完整闭环。
- **自动化交付集成**：**[NEW]** Somatic 模块已正式接入 `DataDeliver` 调度系统。流程会自动解析 `contrasts.csv` 中的配对关系，并自动生成以下核心交付物：
    - `04.filter_variant/somatic/{contrast}.somatic_filtered.vcf.gz` (最终高置信突变集)
    - `04.filter_variant/somatic/{contrast}.somatic_filtered.annotation.csv/html` (SnpEff 功能注释)
    - `03.call_variant/somatic/{contrast}.contamination.table` (交叉污染度评估报告)

### 3. 工程架构优化 (Workflow Management)
- **DataDeliver 调度系统**：采用“模块化定义 + 动态拼接”模式，彻底解耦文件路径定义（`datadeliver.py`）与逻辑调度（`common.py`）。
- **自动化路径映射**：`id_convert.py` 现在支持自动构建 BWA-MEM 规范路径，并精准匹配瘤旁配对样本。
- **规范化监控**：所有规则均对齐了 `06.mapping.smk` 标准格式，包含 `message`、`benchmark` 及 `resources` 资源隔离。

## 🛠️ 核心架构

```text
DNAFlow/
├── rules/
│   ├── 06.mapping.smk       # 比对与质控 (BWA-MEM, MarkDup, Qualimap)
│   ├── 07.callvariant.smk   # Germline 调度器 (支持 Bcf/GATK 动态加载)
│   ├── 07.somatic.smk       # Somatic 流程 (Mutect2 完整过滤链)
│   └── sub_rules/           # 具体算法实现模块 (VQSR, HardFilter 等)
├── utils/
│   ├── common.py            # DataDeliver 总调度逻辑 (分析模式切换)
│   ├── datadeliver.py       # 模块化文件定义库 (各模块产物定义)
│   └── id_convert.py        # 样本表与对比表解析引擎
└── config/
    └── config.yaml          # 流程控制核心配置文件
```

## 📋 快速开始

### 1. 配置分析模式
在 `config.yaml` 中指定你的分析目标：

```yaml
# 分析类型: germline 或 somatic
analysis_type: "somatic"

# 若为 somatic，需指定是否为配对模式
whether_paired_group: True
```

### 2. 准备配对信息 (仅限 Somatic)
创建 `contrasts.csv` 定义 Tumor/Normal 关系（表头需为 Tumor, Normal）：
```csv
Tumor,Normal
Patient01_Tumor,Patient01_Normal
Patient02_Tumor,Patient02_Normal
```

### 3. 运行流程
```bash
# 推荐先进行干跑测试 (Dry run)
snakemake -s snakefile -n

# 正式运行
snakemake -s snakefile --use-conda --cores 32
```

## 🌟 下个版本规划 (Roadmap)
- [ ] 增加 `DeepVariant` AI 调用模块。
- [ ] 引入 `Manta` / `Lumpy` 结构变异 (SV) 检测。
- [ ] 完善基于 Quarto 的变异检测可视化报告。

---

**DNAFlow** - *让每一组碱基序列都清晰可见。*
