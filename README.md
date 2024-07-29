此为 实验物理的大数据方法（2） 课程第二阶段大作业：JUNO probe。

在[一阶段大作业](https://git.tsinghua.edu.cn/physics-data/projects/tpl_junosap)中，我们完成了JUNO探测器的模拟和分析，对JUNO的物理过程有了初步了解，并且通过模拟结果绘制了probe函数的图像。如果你忘了什么是probe函数，请及时回顾，务必明确$r$和$\theta$的定义。

实际上我们第一阶段得到的probe函数只是不含时间的probe函数 $\lambda(r,\theta)$（时间维度被积分了），真正的probe函数是需要考虑时间的，我们通常将含时的probe函数记为$R(r,\theta,t)$，$t$为PE时间。这两个函数之间的关系是：$\lambda(r,\theta)$ = $\int_0^TR(r,\theta,t)\mathrm{d}t$，这里$T$是一个足够大的时间窗口（如600ns或1000ns）。

第二阶段的任务是这样的，我们已经有了足够好的模拟数据，运行 `make data` 即可下载到 data 文件夹中（可以考虑多进程下载，共20个）。你需要用这些训练集，得到一个含时的probe函数$R(r,\theta,t)$。

每个训练集包含10000个顶点，其格式为：

`ParticleTruth` 表
| 名称      | 说明         | 类型    |
| --------- | ------------ | ------- |
| `EventID` | 事件编号          | `'<i4'` |
| `x`       | 顶点坐标x/mm      | `'<f4'` |
| `y`       | 顶点坐标y/mm      | `'<f4'` |
| `z`       | 顶点坐标z/mm      | `'<f4'` |
| `Ek`      | 顶点动能/MeV      | `'<f4'` |
| `Evis`    | 顶点可见能量/MeV  | `'<f4'` |

注：顶点可见能量是指顶点沉积并且用来发光的能量。

`PETruth` 表
| 名称        | 说明           | 类型    |
| ----------- | -------------- | ------- |
| `EventID`   | 事件编号       | `'<i4'` |
| `ChannelID` | PMT 编号       | `'<i4'` |
| `PETime`    | PE 击中时间/ns | `'<f8'` |

几何文件沿用一阶段的`geo.h5`，其格式为：

`Geometry` 表
| 名称        | 说明                     | 类型    |
| ----------- | ------------------------ | ------- |
| `ChannelID` | PMT 编号                 | `'<u4'` |
| `theta`     | 球坐标 $\theta$ **角度** | `'<f4'` |
| `phi`       | 球坐标 $\phi$ **角度**   | `'<f4'` |

也可以通过运行 `make geo.h5` 直接下载。

## 任务说明
本大作业的任务是修改 `probe.py`，造出一个 probe 函数，尽量使分数更高。可以添加额外的文件，但是请不要更改 `coefficient.py` 与 `draw.py` 来作弊，否则判零分。

`probe.py` 中的 `class Probe` 继承了 `coefficient.py` 中的 `class ProbeBase`，因此除了两个抽象函数，其它函数也可以重写。因此你的确可以通过直接重写 `validate` 函数给自己一个高分，这样的情况直接判零分。

## 评分
你可以在本地下载专用的测试文件 `make concat.h5`来自行评分。

画图与评分均使用 `draw.py`。如果需要画图：
```
make draw.pdf
```
如果需要评分：
```
make score
```
这实际上是个 log likelihood，因此可能有负值。分数越大越好。

另外，*为了防止大家对着测试集过拟合，本地评测分数及CI排行仅供参考，最终的黑盒排名将由隐藏测例来决定！*

## 附：评分说明
给定顶点和PMT的相对坐标$(r,\theta)$，我们可以得到一个PE时间序列$\vec{z}=(t_1,t_2,...t_k)$，这里$k$是序列长度。这也是一个非齐次泊松过程，其似然函数为：

$P(\vec{z}|r,\theta)=e^{-\int_0^TR(r,\theta,t)\mathrm{d}t}\prod_{j=1}^kR(r,\theta,t_j)$

如果我们有$N$次采样，那么似然函数就是：

$\mathcal{L}=\prod_{i=1}^NP(\vec{z_i}|r_i,\theta_i)=\prod_{i=1}^Ne^{-\int_0^TR(r_i,\theta_i,t)\mathrm{d}t}\prod_{j=1}^{k_i}R(r_i,\theta_i,t_{ij})$  

取对数得到：

$\log \mathcal{L}=\sum_{i=1}^N(-\int_0^TR(r_i,\theta_i,t)\mathrm{d}t + \sum_{j=1}^{k_i}\log R(r_i,\theta_i,t_{ij}))$

评分时，会将你的probe函数$R(r,\theta,t)$输入到$\log \mathcal{L}$中，得到的值越大越好。

