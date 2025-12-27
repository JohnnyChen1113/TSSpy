# TSSpy 项目开发总结

> 最后更新：2025-12-26
> 版本：v0.11.0
> 匹配率：99.83%

---

## 1. 项目背景

**目标**：用 Python (pysam) 重新实现 R 包 TSSr 的 TSS 提取功能，实现与 TSSr 输出高度一致。

**当前状态**：达到 **99.83%** 匹配率（162,922 / 163,203 位置完全一致）

**GitHub 仓库**：https://github.com/JohnnyChen1113/TSSpy

---

## 2. 关键发现

### 2.1 负链 G-mismatch 移除算法

TSSr 的核心逻辑（之前理解错误，已修复）：

| 步骤 | TSSr 行为 | 说明 |
|------|----------|------|
| 1 | `resize(width=1, fix='start')` | 对负链返回 **END** 位置（不是 start） |
| 2 | `getSeq(Genome, ...)` | 返回参考序列的 **互补碱基** |
| 3 | 检查 `!= 'G'` | 如果互补碱基不是 G，则是 mismatch |

**TSSr 源码位置**：`TSSr/R/ImportFunctions.R` 第 147-157 行

```r
# TSSr 负链 G-mismatch 移除
Gm <- which(substr(seq, start = read.length, stop = read.length) == "C")
while(length(Gm) > 0){
    G.mismatch <- Gm[getSeq(Genome, resize(readsGR.m[Gm], width=1, fix="start")) != "G"]
    end(readsGR.m[G.mismatch]) <- end(readsGR.m[G.mismatch]) - 1
    i = i + 1
    Gm <- G.mismatch[which(substr(seq[G.mismatch], start=1, stop=i) == paste(rep("C",i), collapse=""))]
}
```

### 2.2 序列存储方式

- BAM 文件存储负链 reads 时使用 **正向方向**（不是反向互补）
- Rsamtools 和 pysam 返回 **完全相同** 的序列
- 验证脚本：`test_direction.R`, `test_rsamtools.R`

### 2.3 坐标系统

| 工具 | 坐标系统 | reference_end 含义 |
|------|---------|-------------------|
| pysam | 0-based, half-open | 0-based exclusive |
| TSSr/GRanges | 1-based, closed | 1-based inclusive |

**换算关系**：
- `pysam.reference_start` + 1 = TSSr 的 `start` (1-based)
- `pysam.reference_end` (0-based exclusive) = TSSr 的 `end` (1-based inclusive)

### 2.4 剩余 0.17% 差异来源

TSSr 计算 `mapped.length` 时使用 **CIGAR 数字求和**（包括插入），而非标准的 reference length：

```r
# TSSr 的计算方式 (ImportFunctions.R 第 68 行)
mapped.length <- sum(str_extract_all(cigar, "([0-9]+)"))
# 例如: "39M1I35M" → 39 + 1 + 35 = 75 (而非 reference_length = 74)
```

**影响**：
- 对于有插入的 reads，TSSr 的 end 会比 pysam 的 `reference_end` 大 1
- 差异位置数：TSSpy-only 113 个，TSSr-only 281 个
- 全部差异都在负链

---

## 3. 代码修改历史

### 3.1 迭代记录

| 迭代 | 匹配率 | 主要变化 |
|------|--------|---------|
| 1 | 94.4% | 初始版本，负链算法有误 |
| 2 | 87.5% | 错误修复尝试 |
| 3 | 87.3% | 互补逻辑错误 |
| 4 | 99.83% | 修复 END 位置 + 互补逻辑 |
| 5 | 99.60% | 尝试 query_length（反而更差） |
| 6-7 | 83.6% | 使用错误的参考文件 |
| 7b | 99.83% | 使用正确的参考文件 |

### 3.2 tss_calling.py 核心修改

**版本变化**：0.2.0 → 0.11.0

```python
# 负链 TSS 位置计算（修复后）
if read.is_reverse:
    read_len = len(read_seq)
    pos = read.reference_end  # 使用 pysam 的 reference_end

    # Round 1: 检查最后一个碱基是否是 C
    last_base = read_seq[read_len - 1]
    if last_base != 'C':
        return pos

    # 检查参考基因组在 END 位置的碱基
    ref_end_pos = pos - 1  # 0-based
    ref_base = fasta.fetch(chrom, ref_end_pos, ref_end_pos + 1).upper()

    # TSSr 的 getSeq() 对负链返回互补碱基
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    ref_base_complement = complement.get(ref_base, 'N')

    if ref_base_complement == 'G':
        # 互补碱基是 G（即参考是 C），不是 mismatch
        return pos

    # Mismatch! 减少 TSS 位置
    pos -= 1

    # 后续迭代：检查序列开头是否有连续的 C
    # ...
```

### 3.3 新增模块

| 文件 | 功能 |
|------|------|
| `consensus_cluster.py` | 共识聚类 |
| `merge_samples.py` | 样本合并 |
| `shape_cluster.py` | 形状聚类 |
| `tss_extractor.py` | 备用 TSS 提取 |

---

## 4. 重要文件位置

```
/home/junhaochen/TSSpy/
├── TSSpy/                          # Python 包（Git 仓库）
│   ├── tss_calling.py              # 核心 TSS 提取代码 ★
│   ├── main.py                     # CLI 入口
│   ├── gene_assign.py              # 基因注释
│   ├── consensus_cluster.py        # 共识聚类
│   ├── merge_samples.py            # 样本合并
│   ├── shape_cluster.py            # 形状聚类
│   └── bigwig.py                   # BigWig 导出
│
├── TSSr/                           # R 包源码（参考用）
│   └── R/
│       ├── ImportFunctions.R       # TSSr 核心算法 ★
│       ├── ImportMethods.R         # getTSS 方法
│       └── ...
│
├── TSSpy_development_summary.md    # 本文档
├── iteration_log.md                # 迭代测试日志
├── Rsamtools_analysis.md           # Rsamtools 分析文档
├── BSgenome_analysis.md            # BSgenome 分析文档
│
├── test_*.R                        # R 测试脚本
│   ├── test_direction.R            # 验证序列方向
│   ├── test_getseq.R               # 验证 getSeq 行为
│   ├── test_rsamtools.R            # 验证 Rsamtools 行为
│   ├── test_specific.R             # 特定位置调试
│   └── test_cigar.R                # CIGAR 计算验证
│
├── S01-S04.sorted.bam              # 测试 BAM 文件
├── Scer.genome.fasta               # 参考基因组（符号链接）
├── ALL.samples.TSS.raw.txt         # TSSr 输出（对照用）
└── test_iteration_*.TSS.tsv        # 各迭代测试输出
```

---

## 5. 测试命令

### 5.1 运行 TSSpy

```bash
conda run -n daily python -m TSSpy.main tssCalling \
    -i S01.sorted.bam -i S02.sorted.bam -i S03.sorted.bam -i S04.sorted.bam \
    -o output.TSS.tsv \
    -n "S01 S02 S03 S04" \
    -r Scer.genome.fasta \
    -p 4
```

**注意**：参考基因组文件是 `Scer.genome.fasta`，不是 `sacCer3.fa`

### 5.2 比较结果

```bash
# 排序并比较
cut -f1,2,3 output.TSS.tsv | tail -n +2 | sort > /tmp/tsspy.sorted
cut -f1,2,3 ALL.samples.TSS.raw.txt | tail -n +2 | sort > /tmp/tssr.sorted

# 统计
echo "TSSpy-only:" && comm -23 /tmp/tsspy.sorted /tmp/tssr.sorted | wc -l
echo "TSSr-only:" && comm -13 /tmp/tsspy.sorted /tmp/tssr.sorted | wc -l
echo "Common:" && comm -12 /tmp/tsspy.sorted /tmp/tssr.sorted | wc -l
```

### 5.3 调试特定位置

```python
import pysam

bam = pysam.AlignmentFile('S01.sorted.bam', 'rb')
fasta = pysam.FastaFile('Scer.genome.fasta')

for read in bam.fetch('chrI', 100600, 100700):
    if read.is_reverse:
        print(f'ref_start={read.reference_start}, ref_end={read.reference_end}')
        print(f'CIGAR={read.cigarstring}, seq[-1]={read.query_sequence[-1]}')
```

---

## 6. 待解决问题（如需 100% 匹配）

要完全匹配 TSSr，需要：

1. 对有 **插入 (I)** 的 reads，使用 CIGAR 数字求和计算 end 位置
2. 实现方式：
   ```python
   import re
   def calculate_tssr_mapped_length(cigar_string):
       numbers = re.findall(r'(\d+)', cigar_string)
       return sum(int(n) for n in numbers)

   # 使用 CIGAR sum 代替 reference_end
   tssr_end = read.reference_start + calculate_tssr_mapped_length(read.cigarstring)
   ```

3. **但这在生物学上可能不正确**：插入不消耗参考序列，所以 pysam 的 `reference_end` 更准确

**建议**：当前 99.83% 的匹配率已经足够好，差异仅来自少量有插入的 reads（约 0.6% 的 reads 有 indels）。

---

## 7. Git 信息

- **仓库**：https://github.com/JohnnyChen1113/TSSpy
- **当前版本**：v0.11.0
- **分支**：main

### 提交历史

```bash
git log --oneline -5
# 4302794 v0.11.0: Major update with fixed G-mismatch removal algorithm
# 886411f update the usage of tsspy
# 3a47750 First version of TSSpy
# 35b736b Initial commit
```

---

## 8. 环境信息

```bash
conda activate daily
# Python 3.10
# pysam, pandas, typer
```

---

## 9. 联系与协作

下次继续开发时，可以：
1. 阅读本文档回顾进度
2. 查看 `iteration_log.md` 了解详细迭代过程
3. 使用 `test_*.R` 脚本验证 R 端行为
4. 运行测试命令验证当前代码状态
