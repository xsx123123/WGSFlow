# WGSFlow - Reference Genome Build Pipeline

## 概述

`make_reference` 是 WGSFlow 的参考基因组构建模块，用于自动化生成与 WGSFlow 分析流程完全兼容的参考基因组文件。

## 功能特性

- 自动构建 BWA-MEM 索引
- 生成 samtools faidx 索引
- 创建 GATK 序列字典 (.dict)
- 标准化注释文件 (GTF/GFF/BED12)
- 支持多基因组并行构建
- 与 WGSFlow `config/reference.yaml` 格式完全兼容

## 目录结构

```
make_reference/
├── Snakefile          # 主 Snakemake 流程文件
├── config.yaml        # 配置文件模板
└── README.md          # 本说明文档
```

## 快速开始

### 1. 准备配置文件

复制 `config.yaml` 并根据你的需求修改：

```bash
cd make_reference
cp config.yaml my_config.yaml
# 编辑 my_config.yaml
```

### 2. 配置参考基因组

在配置文件中定义要构建的参考基因组：

```yaml
genomes:
  hg38_p14:
    species: "Homo_sapiens"
    genome_version: "hg38_p14"
    ploidy: 2
    
    input:
      fasta: "/path/to/Homo_sapiens_assembly38.fasta"
      gtf: "/path/to/Homo_sapiens_assembly38.refseq.gtf"
      # ... 其他文件
```

### 3. 运行构建流程

```bash
# 干跑测试
snakemake -s Snakefile --configfile my_config.yaml -n

# 正式运行
snakemake -s Snakefile --configfile my_config.yaml --use-conda --cores 32
```

## 生成的文件结构

构建完成后，会生成以下目录结构：

```
{output_root}/
├── {species}/
│   ├── {genome_version}/
│   │   ├── {genome_version}.fasta          # 基因组序列
│   │   ├── {genome_version}.fasta.fai      # samtools索引
│   │   ├── {genome_version}.dict           # GATK序列字典
│   │   ├── {genome_version}.gtf            # 基因注释
│   │   ├── {genome_version}.gff            # 基因注释
│   │   └── {genome_version}.bed12          # BED12注释
│   └── bwa-mem_v11/
│       ├── {genome_version}.fasta -> ../{genome_version}/{genome_version}.fasta
│       ├── {genome_version}.amb
│       ├── {genome_version}.ann
│       ├── {genome_version}.bwt
│       ├── {genome_version}.pac
│       └── {genome_version}.sa
```

## 配置 WGSFlow 使用新构建的参考基因组

构建完成后，更新 WGSFlow 的 `config/reference.yaml`：

```yaml
reference_path: "/path/to/your/output_root"

can_use_genome_version:
  - hg38_p14
  - your_new_genome

bwa_mem:
  your_new_genome:
    index: {species}/bwa-mem_v11/{genome_version}
    genome_fa: {species}/{genome_version}/{genome_version}.fasta
    genome_gtf: {species}/{genome_version}/{genome_version}.gtf
    # ... 其他配置
```

## 构建选项

在 `config.yaml` 的 `build_options` 部分可以配置：

```yaml
build_options:
  threads: 16                    # 并行线程数
  create_dict: true              # 创建GATK序列字典
  build_bwa: true                # 构建BWA索引
  build_faidx: true              # 构建samtools faidx索引
  generate_bed12: true           # 从GTF生成BED12
```

## 注意事项

1. **输入文件路径**：确保配置文件中的所有输入文件路径都正确且可访问
2. **磁盘空间**：参考基因组构建需要大量磁盘空间（人类基因组约需50GB）
3. **计算资源**：BWA索引构建需要较多内存和CPU时间
4. **物种命名**：建议使用标准的物种拉丁名格式（如 Homo_sapiens）

## 故障排除

### 问题：Conda环境创建失败
确保 WGSFlow 的 `envs/` 目录中有所需的环境文件（bwa.yaml, samtools.yaml, gatk4.yaml）。

### 问题：权限错误
确保输出目录有写入权限，或者修改 `output_root` 为你有权限的目录。

### 问题：内存不足
减少 `build_options.threads` 或在有更多内存的机器上运行。
