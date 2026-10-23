# 独立每细胞供钙当量

## 定义与证据

`C_eq = m f / (M_CaCO3 V)`。一摩尔CaCO3含一摩尔Ca。
M=100.087 g/mol，来源：[NIST Chemistry WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C471341&Mask=69&Units=CAL)
页面顶部 Formula / Molecular weight，CAS 471-34-1（审查时核对）。
独立登记表为 data/processed/calcium_equivalent_parameters.csv；不改主模型108项登记。

pg和pL均带1e-12换算因子，相消后得到mol/L；再乘1e6得µM：

`C_eq[µM] = 1e6 m[pg] f / (100.087 V[pL])`

反算：`m[pg] = C_eq[µM] 100.087 V[pL] / (1e6 f)`，仅在f>0时唯一。
f=0而目标>0不可达；f=0且目标=0对应任意质量，不能返回唯一质量。
正向m=0或f=0结果为0。拒绝负数、非有限数、f>1和V≤0。

## 独立手算例

用户任务书给出的条件：0.1 pg纯CaCO3、50%溶解、2 pL胞质体积。
已溶解质量5e-14 g，Ca物质的量约4.99565e-16 mol，除以2e-12 L，
得到2.49782689e-4 mol/L，即249.782689 µM（约0.25 mM）。
这三个输入是核算例，不是样品实测或文献典型值。

对m、f、V的对数敏感度在正值处分别为+1、+1、-1；
允许的倍数变化只表示公式比例关系，不构造无来源的生物学参数分布。

## 使用方式

在工程根目录：

```powershell
.\.venv\Scripts\python.exe -m cabomb_model.calcium_equivalent --mass-pg 0.1 --dissolved-fraction 0.5 --volume-pl 2 --output-dir outputs/review_round1/calcium
.\.venv\Scripts\python.exe -m cabomb_model.calcium_equivalent --target-um 249.78268906051736 --dissolved-fraction 0.5 --volume-pl 2 --output-dir outputs/review_round1/calcium_inverse
```

可用 `--input-note` 标明细胞系、时间、测量方法或假设依据；不能因填写说明就自动认定实测。
结果包括数值CSV、输入登记、运行输入JSON及中文说明。反算质量是条件推导，不是实验输入。

## 不回答的问题

供钙当量把溶解钙全部按胞质体积折算，不证明钙已转移到胞质，
不等于真实游离Ca2+，不推断死亡阈值或治疗剂量。不引入无依据的“游离比例”。
不能将A_cargo作为钙质量：cargo是独立归一化账本。

整颗粒质量需乘纯CaCO3质量分数；元素钙质量应使用Ca的摩尔质量重新计量，
本接口只接受纯CaCO3。按质量计算不需要颗粒数。由颗粒数推质量另需尺寸、
形状、密度、成分与孔隙假设，当前工具不进行这种换算。

当前IV模型没有绝对剂量到每细胞摄入量的已验证桥接。
培养皿加样浓度也不是每细胞摄入量：还需培养液体积、细胞数、摄取比例和表面/内化区分。
真实游离钙动力学留待后续获得溶解区室、胞质转移、缓冲、内质网/线粒体与外排数据。
