[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Let the Laser Beam Connect the Dots: Forecasting and Narrating Stock Market Volatility

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

This repository contains supporting materials for the paper [Let the Laser Beam Connect the Dots: Forecasting and Narrating Stock Market Volatility](https://doi.org/10.1287/ijoc.2022.0055) by Zhu Zhang, Jie Yuan, and Amulya Gupta.

The software and data in this repository are a snapshot of the software and data that were used in the research reported in the paper.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.xxxx.xxxx

https://doi.org/10.1287/ijoc.xxxx.xxxx.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{xxxx,
  author =        {Zhu Zhang and Jie Yuan and Amulya Gupta},
  publisher =     {INFORMS Journal on Computing},
  title =         {Let the Laser Beam Connect the Dots:Forecasting and Narrating Stock Market Volatility},
  year =          {xxxx},
  doi =           {10.1287/ijoc.xxxx.xxxx.cd},
  url =           {https://github.com/INFORMSJoC/xxxx.xxxx},
}  
```

## Description
We implement LASER in PyTorch. The implementation of LASER-BEAM is based on [constrained_decoding](https://github.com/chrishokamp/constrained_decoding/tree/master).

This is an example of how to set up the project locally and run [LASER](results/LASER4.png) model.

## Building

1. Clone and install the repo
   ```
   git clone https://github.com/INFORMSJoC/2022.0055.git
   ```
  
2. Create an environment using
   ```
   conda create --name LASER --file requriements.txt
   ```

## Results

[Table 1](results/Table1.png) shows the results for market volatility forecasting. Lower MAPE and higher Precision, Recall, F-1 are better.

[Table 2](results/Table2.png) shows the performance evaluation for model narration. (All scores are average over the collection of narratives.
Higher Fluency, Informativeness, and Harmonic Mean are better)

[Table 3](results/Table3.png) shows the time complexity of each component in LASER.

[Table 6](results/Table6.png) shows the instance-level model performance when H = 1.

[Table 7](results/Table7.png) shows an example: Generated narratives for short-term thread (forecasting model is LASER, and horizon H = 1).

[Table 8](results/Table8.png) shows an example: Generated narratives for long-term threads (forecasting model is LASER model, and horizon H = 1).

[Figure 6](results/F6.png) shows the qualitative analysis of volatility forecasting, short horizon (H=1).

## Replicating

<h3 align="left">Generate datasets</h3>
Download data from https://www.wsj.com/news/archive/.

To replicate the experiment, first, generate datasets using pre-trained BERT

  ```
  python scripts/generate_dataset.py
  ```

<h3 align="left">Run training/validation/testing</h3>
There are two tasks: SPVt1_future and SPVt22_future.

* Run SPVt1_future: 
  ```
  python src/train_val_test.py SPVt1_future
  ```

* Run SPVt22_future:
  ```
  python src/train_val_test.py SPVt22_future
  ```
