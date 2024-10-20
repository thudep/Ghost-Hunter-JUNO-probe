from multiprocessing import Pool
import argparse
import numpy as np
import h5py
from tqdm import tqdm

# PMT数量
N = 17612
# 液闪区域半径(mm)
R0 = 17710
# 数据文件编号
Seq = np.arange(16001, 16002)
# 数据文件个数
N0 = 1
# probe函数积分上限(ns)
T_max = 1000

def get_probe(Arguments):
    """
    get_probe() 函数

    本函数用于对一组数据计算probe函数。

    输入为数据文件序号i，数据文件夹名data，几何数据geo_theta和geo_phi，
    r和theta方向bins格子数量和theta方向格子数量。
    """
    i, argument, geo_theta, geo_phi = Arguments

    # read the file
    with h5py.File(f"{argument[0]}/{i}.h5", 'r') as h5file_r:
        ParticleTruth = h5file_r["ParticleTruth"][...]
        PETruth = h5file_r["PETruth"][...]

    # 计算各顶点的归一化位置, 和各PMT的位置
    x = ParticleTruth['x'][:, None] / R0
    y = ParticleTruth['y'][:, None] / R0
    z = ParticleTruth['z'][:, None] / R0

    n = np.size(ParticleTruth)
    probe = np.zeros((argument[1], argument[1], argument[2]))

    # 计算各[顶点, PMT]对应的probe函数的r和θ
    r = np.sqrt(x**2 + y**2 + z**2) * np.ones((n, N))
    θ = np.arccos(np.clip((x * np.sin(geo_theta) * np.cos(geo_phi)
                            + y * np.sin(geo_theta) * np.sin(geo_phi)
                            + z * np.cos(geo_theta)) / r, -1, 1)).flatten()
    θ = np.min(np.array([θ, 2*np.pi - θ]), axis=0)
    r = r.flatten()

    # 对每个t格子
    for j in tqdm(range(argument[2])):
        Bin_T = np.arange(j, j+2) * T_max / argument[2]

        # 每个(r, θ)组合出现的次数统计
        hist_pmt_j = np.histogramdd(
            np.array([PETruth["EventID"], PETruth["ChannelID"], PETruth["PETime"]]).T,
            bins=[np.arange(len(ParticleTruth)+1), np.arange(N+1), Bin_T])[0].flatten()

        # 在每个格子里对(r, θ)组合出现的次数取平均, 结果为(sum_probe / num)
        sum_probe = np.histogram2d(r, θ,
            bins=[argument[3], argument[4]], weights=hist_pmt_j)[0]
        num = np.histogram2d(r, θ, bins=[argument[3], argument[4]])[0]

        probe[:, :, j] += np.where(num > 0, sum_probe / num, 0)

        # 将数组的零值替换为其中其余元素的平均值
        probe[:, :, j][probe[:, :, j] == 0] = probe[:, :, j][probe[:, :, j] != 0].mean()
    
    return probe
def main():
    """
	main() 函数

	此函数用于计算JUNO探测器的三维probe函数（r, θ, t），并将计算结果保存为HDF5文件。
    输入数据包括几何文件和多个训练数据文件。

	函数功能：
	1. 解析命令行参数，获取输入的几何文件、训练数据文件夹路径、输出文件路径、空间和时间的分箱数。
	2. 读取几何文件，提取PMT的球坐标角度θ和φ。
	3. 依次读取每个训练数据文件，计算各顶点和PMT的相对位置和对应的r和θ值。
	4. 统计每个(r, θ, t)组合出现的次数，并对所有训练文件进行累加。
	5. 将累加的探测器响应函数存储为HDF5格式，并将分箱信息作为属性保存。

	参数：
	- `-g, --geo`: 几何文件路径（HDF5格式），包含PMT的球坐标角度信息。
	- `--data`: 训练数据文件夹路径，包含多个训练数据文件（HDF5格式）。
	- `-o, --output`: 输出文件路径（HDF5格式），用于保存计算得到的探测器响应函数。
	- `-b, --bins`: 空间分箱数，用于对r和θ进行分箱。
	- `-t, --tbins`: 时间分箱数，用于对时间t进行分箱。

	输入：
	- 几何文件（HDF5格式）：包含PMT的球坐标角度信息（θ和φ）。
	- 训练数据文件（HDF5格式）：包含顶点和光电子信息。

	输出：
	- 探测器响应函数的三维数组（r, θ, t），并以HDF5格式保存到指定的输出文件中。
    输出文件中还包含分箱信息作为属性保存。

	注意事项：
	- 确保输入的几何文件和训练数据文件路径正确无误。
	- 输出文件路径应为有效的HDF5文件路径。

	示例用法：
	```bash
	./histogram.py -g geometry.h5 --data data_folder -o output_probe.h5 -b 20 -t 100
	```
    """

    # Argument parsing
    psr = argparse.ArgumentParser()
    psr.add_argument("-g", "--geo", dest="geo", type=str, help="geometry file")
    psr.add_argument("--data", dest="data", type=str, help="training data")
    psr.add_argument("-o", "--output", dest="opt", type=str, help="output file")
    psr.add_argument("-b", "--bins", dest="Bins", type=str, help="output file")
    psr.add_argument("-t", "--tbins", dest="T_Bins", type=str, help="output file")
    args = psr.parse_args()

    # read geometry data
    with h5py.File(args.geo, 'r') as h5file_r:
        geo = h5file_r["Geometry"][...]
    geo_theta = np.deg2rad(geo["theta"][None, :N])
    geo_phi = np.deg2rad(geo["phi"][None, :N])

    r_bins = (np.arange(int(args.Bins)+1) / int(args.Bins)) ** (1/3)
    theta_bins = np.arccos(np.arange(int(args.Bins)+1) / int(args.Bins))[::-1]

    with Pool(processes=20) as pool:
        probe = np.sum(np.array(pool.map(get_probe,
            zip(Seq, zip([args.data]*N0, [int(args.Bins)]*N0,
            [int(args.T_Bins)]*N0, [r_bins]*N0, [theta_bins]*N0),
            [geo_theta]*N0, [geo_phi]*N0))), axis=0)

    probe = probe / N0 * int(args.T_Bins) / T_max

    with h5py.File(args.opt, 'w') as h5file_w:
        dataset = h5file_w.create_dataset('Probe', data=probe)
        dataset.attrs['Bins'] = int(args.Bins)
        dataset.attrs['T_Bins'] = int(args.T_Bins)
        dataset.attrs['R_Bins'] = r_bins
        dataset.attrs['Theta_Bins'] = theta_bins

if __name__ == "__main__":
    main()
