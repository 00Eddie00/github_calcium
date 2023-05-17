import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


def cal_temporal_properties(ampl, time_interval):
    peak_value = np.max(ampl)  # 峰值，即某空间点在时间尺度上的最大荧光强度
    t_rise_index = np.argmax(ampl)  # 峰化时间对应ampl中的索引

    temp = np.maximum(ampl - peak_value / 2, peak_value / 2 - ampl)
    t50_index = np.argmin(temp[t_rise_index:])  # 从最大衰减到50%时索引
    fdhm_index = np.argmin(temp[:t_rise_index])  # 从50%增加到最大时索引

    t_rise = t_rise_index * time_interval  # 峰化时间，即荧光强度从基线达到峰值所需时间
    t50 = t50_index * time_interval  # 荧光强度从峰值衰减50%所需时间
    fdhm = (t50_index + t_rise_index - fdhm_index) * time_interval  # full duration at half maximum，半峰全持续时间
    fdhm_start = fdhm_index * time_interval

    return peak_value, t_rise, t50, fdhm, fdhm_start


def temporal_plotter(temporal_path, time_interval, data_source, save=False):
    temporal_fluorescence = np.loadtxt(temporal_path)

    ampl = temporal_fluorescence / temporal_fluorescence[0] - 1  # 振幅，即荧光强度：fluorescence intensity
    t = [i * time_interval for i in range(len(temporal_fluorescence))]  # 时间，单位ms
    peak_value, t_rise, t50, fdhm, fdhm_start = cal_temporal_properties(ampl, time_interval)  # 计算钙火花的各属性

    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(t, ampl, ls='-', lw=2)

    # 绘制tRise时间点及峰值
    plt.plot(t_rise, peak_value, 'o', color='orange')
    plt.text(t_rise, peak_value, f'tRise={np.round(t_rise, 4)}ms peak={np.round(peak_value, 4)}', fontsize=14)
    # 绘制t50时间点
    t50_t_index = t50 + t_rise
    t50_ampl = ampl[int((t50 + t_rise) / time_interval)]
    plt.plot(t50_t_index, t50_ampl, 'o', color='orange')
    plt.text(t50_t_index, t50_ampl, f't50={np.round(t50, 4)}ms', fontsize=14)
    # 绘制FDHM
    fdhm_start_ampl = ampl[int(fdhm_start / time_interval)]
    plt.plot(fdhm_start, fdhm_start_ampl, 'o', color='orange')
    plt.plot([fdhm_start, t50_t_index], [fdhm_start_ampl, t50_ampl], ls='--', color='orange')
    plt.text(t_rise, fdhm_start_ampl, f'FDHM={np.round(fdhm, 4)}ms', fontsize=12)
    plt.text(0, 0, f'FDHM={np.round(fdhm, 4)}ms', fontsize=12)

    plt.xlabel('t(ms)', fontsize=14)
    plt.ylabel('ΔF/F0', fontsize=14)
    plt.suptitle('中心：时间分布')
    plt.title(f"数据来源：{data_source}")

    plt.grid()
    if save:
        plt.savefig(f"{temporal_path}.png")
    plt.show()


def cal_spatial_properties(spatial_fluorescence, dis_interval):
    radius = (len(spatial_fluorescence) - 1) / 2 * dis_interval

    max_fluo = np.max(spatial_fluorescence)
    min_fluo = np.min(spatial_fluorescence)
    half_maximum = (max_fluo + min_fluo) / 2
    temp = np.maximum(spatial_fluorescence - half_maximum, half_maximum - spatial_fluorescence)

    fwhm_start = np.argmin(temp) * dis_interval - radius
    fwhm_end = -fwhm_start
    fwhm = 2 * fwhm_end  # full width at half maximum，半峰全宽

    return fwhm, fwhm_start, fwhm_end


def spatial_plotter(spatial_path, dis_interval, save=True):
    spatial_fluorescence = np.loadtxt(spatial_path)

    radius_count = (len(spatial_fluorescence) - 1) / 2
    x = [(i - radius_count) * dis_interval for i in range(len(spatial_fluorescence))]

    fwhm, fwhm_start, fwhm_end = cal_spatial_properties(spatial_fluorescence, dis_interval)  # 计算钙火花的fwhm

    plt.figure(dpi=300, figsize=(10, 10))
    plt.plot(x, spatial_fluorescence, ls='-', lw=2)
    # 绘制FWHM
    fwhm_start_fluo = spatial_fluorescence[int(fwhm_start / dis_interval + radius_count)]
    fwhm_end_fluo = spatial_fluorescence[int(fwhm_end / dis_interval + radius_count)]
    plt.plot(fwhm_start, fwhm_start_fluo, 'o', color='orange')
    plt.plot(fwhm_end, fwhm_end_fluo, 'o', color='orange')
    plt.plot([fwhm_start, fwhm_end], [fwhm_start_fluo, fwhm_end_fluo], ls='--', color='orange')
    plt.text(0, fwhm_start_fluo, f'FWHM={np.round(fwhm, 4)}nm', fontsize=12)

    plt.xlabel('x(nm)', fontsize=14)
    plt.ylabel('concentration(mM)', fontsize=14)
    plt.suptitle('空间分布')
    data_source = spatial_path[spatial_path.find("ResultSet_"):]
    plt.title(data_source)

    plt.grid()
    if save:
        plt.savefig(f"{spatial_path}.png")
    plt.show()


if __name__ == '__main__':
    my_temporal_path = "E:\\Study\\sr&sl参数调整\\情况九子情况3\\NANO_CaG_avg_c.csv"
    my_time_interval = 2 * 10 ** -6 * 100 * 1000
    temporal_plotter(my_temporal_path, my_time_interval, "情况九子情况3-CaG浓度时间变化", save=True)

    # spatial_path = "E:\\Study\\Fixed Unstructured ResultSets\\ResultSet_22-0518-1459_只取上次迭代的值\\CaG_psf_spatial_distribution.csv"
    # my_dis_interval = 1
    # spatial_plotter(spatial_path, my_dis_interval)
