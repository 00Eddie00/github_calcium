# sr 、 sl 参数调整

## 参数情况

初始参数：
参考论文：1998 Smith钙火花论文

- sr membrane
    - K_PLUS_SRM = 115000
    - K_MINUS_SRM = 100
    - TOTAL_SRM = 0.047
    - K_DIV_SRM = K_MINUS_SRM / K_PLUS_SRM = 100 / 115000 = 8.6957*10**-4

- sl membrane
    - K_PLUS_SLM = 115000
    - K_MINUS_SLM = 1000
    - TOTAL_SLM = 1.124
    - K_DIV_SLM = K_MINUS_SLM / K_PLUS_SLM = 1000 / 115000 = 8.6957*10**-3

情况一：

- 其它变化：	
  - 修改Kryr（钙空穴的kryr*2）
  - 修改 calmodulin、troponin c （与1998 Peskoff 论文相同）
    参考论文：1998 Peskoff 论文

- sr membrane
    - K_MINUS_SRM = 100
    - TOTAL_SRM = 13
    - K_DIV_SRM = 0.013
    - K_PLUS_SRM = K_MINUS_SRM / K_DIV_SRM = 100 / 0.013 = 7692.307

- sl membrane
    - K_MINUS_SLM = 1000
    - TOTAL_SLM = 165
    - K_DIV_SLM = 1.1
    - K_PLUS_SLM = K_MINUS_SLM / K_DIV_SLM = 1000 / 1.1= 909.0909

    ![test](.\assests\NANO_Ca_avg_c_情况一.png)

    ![test](.\assests\NANO_CaG_avg_c_情况一.png)

情况二：
将情况一的 TOTAL_SRM、TOTAL_SLM 调为十分之一，共三种子情况：

1. 只调 TOTAL_SRM = 1.3
2. 只调 TOTAL_SLM = 16.5
3. 同时调 TOTAL_SRM = 1.3、TOTAL_SLM = 16.5

---

情况三：
在情况二的基础上第 2 条的基础上，将 K_MINUS_SRM 调为 1000

- sr membrane
    - K_MINUS_SRM = 1000
    - TOTAL_SRM = 13
    - K_DIV_SRM = 0.013
    - K_PLUS_SRM = K_MINUS_SRM / K_DIV_SRM = 1000 / 0.013 = 76923.07

- sl membrane
    - K_MINUS_SLM = 1000
    - TOTAL_SLM = 16.5
    - K_DIV_SLM = 1.1
    - K_PLUS_SLM = K_MINUS_SLM / K_DIV_SLM = 1000 / 165 = 909.0909


情况四：
结合Smith钙火花论文与 1998 Peskoff 论文。
K_PLUS、K_MINUS、K_DIV 使用Smith钙火花论文参数，TOTAL 使用 Peskoff 论文参数
- sr membrane
    - K_PLUS_SRM = 115000
    - K_MINUS_SRM = 100
    - TOTAL_SRM = 13
    - K_DIV_SRM = K_MINUS_SRM / K_PLUS_SRM = 100 / 115000 = 8.6957*10**-4
- sl membrane
    - K_PLUS_SLM = 115000
    - K_MINUS_SLM = 1000
    - TOTAL_SLM = 165
    - K_DIV_SLM = K_MINUS_SLM / K_PLUS_SLM = 1000 / 115000 = 8.6957*10**-3

不行，会溢出

情况五：
在四的基础上，TOTAL_SLM = 16.5
- sr membrane
    - K_PLUS_SRM = 115000
    - K_MINUS_SRM = 100
    - TOTAL_SRM = 13
    - K_DIV_SRM = K_MINUS_SRM / K_PLUS_SRM = 100 / 115000 = 8.6957*10**-4

- sl membrane
    - K_PLUS_SLM = 115000
    - K_MINUS_SLM = 1000
    - TOTAL_SLM = 16.5
    - K_DIV_SLM = K_MINUS_SLM / K_PLUS_SLM = 1000 / 115000 = 8.6957*10**-3

不行，有负值

![test](.\assests\NANO_Ca_avg_c_不包含情况四.png)

![test](.\assests\NANO_Ca_avg_c_不包含情况四、五.png)



![test](.\assests\NANO_CaG_avg_c_不包含情况四.png)

![test](.\assests\NANO_CaG_avg_c_不包含情况四、五.png)

情况六：
1. 在情况二-子情况2基础上调整K_OUT，原来的10倍，即7.854*10**6
2. 在情况三的基础上调整K_OUT，原来的10倍，即7.854*10**6

![NANO_Ca_avg_c_情况六_子情况1](.\assests\NANO_Ca_avg_c_情况六_子情况1.png)![NANO_CaG_avg_c_情况六_子情况1](.\assests\NANO_CaG_avg_c_情况六_子情况1.png)

![NANO_Ca_avg_c_情况六_子情况2](.\assests\NANO_Ca_avg_c_情况六_子情况2.png)

![NANO_CaG_avg_c_情况六_子情况2](.\assests\NANO_CaG_avg_c_情况六_子情况2.png)

![NANO_Ca_avg_c_情况六对比](.\assests\NANO_Ca_avg_c_情况六对比.png)

![NANO_CaG_avg_c_情况六对比](.\assests\NANO_CaG_avg_c_情况六对比.png)

情况七：

1. 在情况二-子情况2基础上调整K_OUT，原来的20倍，即2×7.854×10**6
2. 在情况三的基础上调整K_OUT，原来的20倍，即2×7.854×10**6
3. 在情况二-子情况2基础上调整K_OUT，原来的50倍，即5×7.854×10**6
4. 在情况三的基础上调整K_OUT，原来的50倍，即5×7.854×10**6![NANO_Ca_avg_c_情况六与情况七对比](.\assests\NANO_Ca_avg_c_情况六与情况七对比.png)

![NANO_CaG_avg_c_情况六与情况七对比](.\assests\NANO_CaG_avg_c_情况六与情况七对比.png)

情况八：

在 K_OUT 扩大10倍、KRYR为2 * 6.5 * 10 ** 7的基础上，分别调整 TOTAL_SRM、K_DIV_SRM、TOTAL_SLM、K_DIV_SLM 变小。

- 情况一：TOTAL_SRM=0.13，即缩小至百分之一
- 情况二：K_DIV_SRM=0.0013，即所缩小至十分之一
- 情况三：TOTAL_SLM=1.65，即缩小至百分之一
- 情况四：K_DIV_SLM=0.11，即缩小至十分之一

![NANO_Ca_avg_c_情况二与情况八对比](.\assests\NANO_Ca_avg_c_情况二与情况八对比.png)

![NANO_CaG_avg_c_情况二与情况八对比](.\assests\NANO_CaG_avg_c_情况二与情况八对比.png)

情况九：在情况八的子情况1的基础上调整荧光染料CaG的参数

- 子情况1：K_PLUS_CAG 缩小至 1/5，即 27000 / 5 = 5400

- 子情况2：K_MINUS_CAG 扩大至 5 倍即17 * 5 = 85

- 子情况3：TOTAL_CAG 扩大至 10 倍，即 0.01 * 10 = 0.1


钙离子的变化不大，有重叠曲线。

![NANO_Ca_avg_c_情况九对比](.\assests\NANO_Ca_avg_c_情况九对比.png)

情况九子情况2的T50有明显降低，接近50ms。

![NANO_CaG_avg_c_情况九对比](.\assests\NANO_CaG_avg_c_情况九对比.png)

![NANO_CaG_avg_c_情况九子情况1与2对比](.\assests\NANO_CaG_avg_c_情况九子情况1与2对比.png)

情况十：网格直径450nm，D_CA = 7.0 * 10 ** 8，K+=50000，K-=34

- 子情况2：情况九子情况3
