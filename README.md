# 2024 Ghost Hunter 的赛题: JUNO probe

在宣讲会上, 我们已经对 JUNO 的物理过程有了初步了解, 并且了解了 probe 函数的定义. probe 函数用来描述探测器对顶点的响应, 我们通常将含时的 probe 函数记为 $R(r,\theta,t)$, 该函数表示在相对坐标 $(r, \theta)$ 下随时间 $t$ 变化的 PMT 上接收到的期望 PE 数. 将时间维度积分后, 得到不含时间的 probe 函数 $\lambda(r,\theta)$. 即: 这两个函数之间的关系是: $\lambda(r,\theta)$ = $\int_0^TR(r,\theta,t)\mathrm{d}t$, 这里 $T$ 是一个足够大的时间窗口, 我们取 1000ns.

赛题是这样的, 我们已经有了足够好的模拟数据作为训练集, 你需要用这些训练集, 得到一个 probe 函数.

每个训练集包含10000个顶点, 其格式为:

h5文件中的`ParticleTruth` 表
| 名称      | 说明             |
| --------- | ---------------- |
| `EventID` | 事件编号         |
| `x`       | 顶点坐标x/mm     |
| `y`       | 顶点坐标y/mm     |
| `z`       | 顶点坐标z/mm     |
| `Ek`      | 顶点动能/MeV     |
| `Evis`    | 顶点可见能量/MeV |

注: 顶点可见能量是指顶点沉积并且用来发光的能量.

parquet文件中的
| 名称        | 说明                |
| ----------- | ------------------- |
| `evtID`     | 事件编号            |
| `pmtID`     | PMT 编号            |
| `hitTime`   | PE 击中时间/ns      |
| `lightTime` | 光子被发出的时间/ns |

几何文件`geo.h5`, 其格式为:

`Geometry` 表
| 名称        | 说明                     |
| ----------- | ------------------------ |
| `ChannelID` | PMT 编号                 |
| `theta`     | 球坐标 $\theta$ **角度** |
| `phi`       | 球坐标 $\phi$ **角度**   |

## 任务说明

任务是修改 `probe.py`, 造出一个 probe 函数, 尽量使分数更高. 可以添加额外的文件, 但是请不要更改 `coefficient.py` 与 `draw.py` 来作弊, 否则成绩无效.

`probe.py` 中的 `class Probe` 继承了 `coefficient.py` 中的 `class ProbeBase`, 因此除了两个抽象函数, 其它函数也可以重写. 因此你的确可以通过直接重写 `validate` 函数给自己一个高分, 这样的情况直接判零分.

我们将提供一个测试集供大家在本地评分作为参考. 但是, **为了防止大家对着测试集过拟合, 本地评测分数及CI排行仅供参考, 最终的排名将由隐藏测例来决定!**

## 附: 评分说明

给定顶点和 PMT 的相对坐标 $(r,\theta)$, 我们可以得到一个 PE 时间序列 $\vec{z}=(t_1,t_2,...t_k)$, 这里 $k$ 是序列长度. 这也是一个非齐次泊松过程, 其似然函数为:

$P(\vec{z}|r,\theta)=e^{-\int_0^TR(r,\theta,t)\mathrm{d}t}\prod_{j=1}^kR(r,\theta,t_j)$

如果我们有 $N$ 次采样, 那么似然函数就是:

$\mathcal{L}=\prod_{i=1}^NP(\vec{z_i}|r_i,\theta_i)=\prod_{i=1}^Ne^{-\int_0^TR(r_i,\theta_i,t)\mathrm{d}t}\prod_{j=1}^{k_i}R(r_i,\theta_i,t_{ij})$

取对数得到:

$\log \mathcal{L}=\sum_{i=1}^N(-\int_0^TR(r_i,\theta_i,t)\mathrm{d}t + \sum_{j=1}^{k_i}\log R(r_i,\theta_i,t_{ij}))$

评分时, 会将你的 probe 函数 $R(r,\theta,t)$ 输入到 $\log \mathcal{L}$ 中, 得到的值越大越好.
