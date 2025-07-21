import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_histogram_from_csv(
    csv_file_path,
    column_name,
    num_bins=None,
    title=None,
    xlabel=None,
    ylabel="Frequency", # 默认Y轴标签改回Frequency或Count
    xlim_min=None,
    xlim_max=None,
    ylim_min=None,
    ylim_max=None,
    bar_color='skyblue', # 新增：直方图柱子的颜色
    show_mean=False, # 新增：是否显示平均值
    save_path=None
):
    """
    读取CSV文件的指定列，绘制其频率分布直方图，并可选地显示平均值和自定义颜色。

    Args:
        csv_file_path (str): CSV 文件的完整路径。
        column_name (str): 要绘制直方图的列的名称。
        num_bins (int, optional): 直方图的 bin 数量。
        title (str, optional): 直方图的标题。
        xlabel (str, optional): X 轴的标签。
        ylabel (str, optional): Y 轴的标签。默认为 "Frequency"。
        xlim_min (float/int, optional): X 轴的最小值。
        xlim_max (float/int, optional): X 轴的最大值。
        ylim_min (float/int, optional): Y 轴的最小值。
        ylim_max (float/int, optional): Y 轴的最大值。
        bar_color (str, optional): 直方图柱子的颜色。默认为 'skyblue'。
        show_mean (bool, optional): 是否在图上显示平均值。默认为 False。
        save_path (str, optional): 保存图表的完整路径。
    """
    print("--- 开始绘制直方图 ---")
    print("CSV 文件: {}".format(csv_file_path))
    print("目标列: '{}'".format(column_name))
    if show_mean:
        print("将显示平均值。")

    if not os.path.exists(csv_file_path):
        print("错误: CSV 文件 '{}' 不存在。请检查路径。".format(csv_file_path))
        return

    try:
        df = pd.read_csv(csv_file_path)
    except Exception as e:
        print("错误: 读取 CSV 文件失败：{}".format(e))
        return

    if column_name not in df.columns:
        print("错误: CSV 文件中未找到列 '{}'。".format(column_name))
        print("可用列: {}".format(df.columns.tolist()))
        return

    data = pd.to_numeric(df[column_name], errors='coerce')
    data = data.dropna()

    if data.empty:
        print("警告: 目标列 '{}' 中没有有效的数值数据可用于绘制直方图。".format(column_name))
        return

    print("找到 {} 个有效数值数据点。".format(len(data)))

    plt.figure(figsize=(10, 6))
    
    plt.hist(data, bins=num_bins, edgecolor='black', alpha=0.7, color=bar_color) 

    if title is None:
        plt.title("Frequency Distribution of '{}'".format(column_name))
    else:
        plt.title(title)
        
    if xlabel is None:
        plt.xlabel(column_name)
    else:
        plt.xlabel(xlabel)
        
    plt.ylabel(ylabel)

    plt.grid(axis='y', alpha=0.75)

    if xlim_min is not None and xlim_max is not None:
        plt.xlim(xlim_min, xlim_max)
        print("X轴范围已设置为: [{}, {}]".format(xlim_min, xlim_max))
    
    if ylim_min is not None and ylim_max is not None:
        plt.ylim(ylim_min, ylim_max)
        print("Y轴范围已设置为: [{}, {}]".format(ylim_min, ylim_max))

    if show_mean:
        mean_value = data.mean()
        _, max_y = plt.ylim() 
        
        plt.axvline(mean_value, color='red', linestyle='dashed', linewidth=2, label='Mean')
        
        plt.text(mean_value + (plt.xlim()[1] - plt.xlim()[0]) * 0.01, # x轴范围的1%作为偏移
                 max_y * 0.8, # Y轴80%高度
                 'Mean: {:.2f}'.format(mean_value), 
                 color='red', 
                 horizontalalignment='left',
                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
        plt.legend()
    if save_path:
        try:
            plt.savefig(save_path)
            print("直方图已保存至: '{}'".format(save_path))
        except Exception as e:
            print("错误: 保存图表到 '{}' 失败：{}".format(save_path, e))
        plt.close()
    else:
        print("直方图已生成并显示。")
        plt.show()

    print("--- 直方图绘制完成 ---")


if __name__ == "__main__":
    # --- 用户配置区 ---
    # 你的 CSV 文件路径
    CSV_FILE_PATH = "/mnt/data/renxinyu/rosetta/Tools/3rd_proteinMPNN.csv" 
    SAVE_PLOT_AS = "3rd_dg_proteinMPNN.png" 
    NUMBER_OF_BINS = 20
    Y_AXIS_LABEL = "Frequency" 
    SHOW_MEAN_VALUE = True # 设置为 True 来显示平均值


    # 必须修改的/////////////////////////////////////////////////////////////////////
    TARGET_COLUMN_NAME = "dG_cross" 
    # TARGET_COLUMN_NAME = "total_score" 

    PLOT_TITLE = "dG of proteinMPNN" 
    # PLOT_TITLE = "dG of ProteinMPNN" 
    # PLOT_TITLE = "Total Score of HighMPNN" 

    # X_AXIS_LABEL = "Total Score" 
    X_AXIS_LABEL = "dG"

    # HISTOGRAM_BAR_COLOR = '#e7bfc0' # red
    # HISTOGRAM_BAR_COLOR = '#5b66a1' # blue
    HISTOGRAM_BAR_COLOR = '#bcbcbc' # gray
    

    X_LIMIT_MIN = None
    X_LIMIT_MAX = None 
    Y_LIMIT_MIN = None
    Y_LIMIT_MAX = None
    

    # --- 调用主函数执行任务 ---
    plot_histogram_from_csv(
        csv_file_path=CSV_FILE_PATH,
        column_name=TARGET_COLUMN_NAME,
        num_bins=NUMBER_OF_BINS,
        title=PLOT_TITLE,
        xlabel=X_AXIS_LABEL,
        ylabel=Y_AXIS_LABEL,
        xlim_min=X_LIMIT_MIN,
        xlim_max=X_LIMIT_MAX,
        ylim_min=Y_LIMIT_MIN,
        ylim_max=Y_LIMIT_MAX,
        bar_color=HISTOGRAM_BAR_COLOR,
        show_mean=SHOW_MEAN_VALUE,
        save_path=SAVE_PLOT_AS
    )

    print("\n脚本执行完毕。")