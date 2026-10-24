# CABOMB 个人建模项目

作者：高书伦。包含为 BNU-China 2026 iGEM CABOMB 项目开发的建模代码，以及后续个人研究扩展。

原模型用常微分方程描述 CaCO₃ 载体从血液循环到肿瘤入胞、降解和释放的过程，比较未修饰、OCN、YSA 和双修饰四种设计。后续分析研究早期释放曲线能否预测后期，以及入胞后的等待时间怎样影响累计响应。

## 交互可视化

![原交互展示数据的曲线回放](visualizations/organ_replay/time_replay.gif)

[三维时间回放及源码](visualizations/organ_replay/README.md)保留了之前做的肝脏、肿瘤组织和细胞内释放展示。下载其中的 `index.html` 后，用浏览器打开即可播放、拖动时间和旋转模型。这一页使用原来的 V1.2、100 nm OCN–YSA 情景；GIF 是本次用原展示的保存数据导出的预览。肝脏参考网格来自 HuBMAP，素材署名和许可证一并保留。

## Usage

在本目录打开 PowerShell，使用 Python 3.10 或更新版本：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install ".[dev,app]"
python -m cabomb_model validate-registry
python -m cabomb_model run-baseline --model-version v1.3 --output-dir outputs/baseline
python -m cabomb_release synthetic --output-dir outputs/synthetic
python -m cabomb_release bridge --output-dir outputs/bridge
python -m pytest
```

如果不需要界面，可以把安装命令中的 `[dev,app]` 换成 `[dev]`。`synthetic` 和 `bridge` 需要使用尚不存在的输出目录，重复运行时换一个目录名即可。

输入参数在 `data/processed/`，程序生成的数据和图在指定的 `outputs/` 子目录里。README 只说明运行方法，不放输入或结果数据表。文献原始附件不随包发布，来源和自备数据的运行方式见 [数据来源](DATA_SOURCES.md) 和 [运行说明](docs/run_guide.md)。

目前这个版本适合展示模型实现、设计比较和方法分析。已有结果及后续值得补的实验见 [研究记录](research/findings.md)。

## AI 辅助说明

本项目由高书伦开发和维护，最初的 ODE 框架和核心代码由本人从零搭建。后续方法探索和项目整理过程中，使用了 AI 工具辅助文献检索、方案讨论、代码编写、调试和文档整理，也使用另一 AI 代理检查部分代码和结果。

AI 生成的内容可能存在错误，研究结论应结合原始文献、代码和验证记录判断。AI 辅助检查不等同于人工同行评审。

项目保留 CABOMB 团队背景说明；引用的文献及公开数据归原作者所有，并单独注明来源。
