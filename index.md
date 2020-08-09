## Welcome to CacheTest

This is the documentation for CacheTest

## Building

In Linux, to build the version that multiplies all elements of a vector by a
constant (used to obtain the results in [Figure 1](https://github.com/INFORMSJoC/JoCTemplate/blob/master/results/mult-test.png) in the
paper), stepping K elements at a time, execute the following commands.

```
make mult
```

Alternatively, to build the version that sums the elements of a vector (used
to obtain the results [Figure 2](https://github.com/INFORMSJoC/JoCTemplate/blob/master/results/sum-test.png) in the paper), stepping K
elements at a time, do the following.

```
make clean
make sum
```

Be sure to make clean before building a different version of the code.

## Results

Figure 1 in the paper shows the results of the multiplication test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](https://github.com/INFORMSJoC/JoCTemplate/blob/master/results/mult-test.png)

Figure 2 in the paper shows the results of the sum test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](https://github.com/INFORMSJoC/JoCTemplate/blob/master/results/sum-test.png)

## Replicating

To replicate the results in [Figure 1](https://github.com/INFORMSJoC/JoCTemplate/blob/master/results/mult-test.png), do either

```
make mult-test
```
or
```
python test.py mult
```
To replicate the results in [Figure 2](https://github.com/INFORMSJoC/JoCTemplate/blob/master/results/sum-test.png), do either

```
make sum-test
```
or
```
python test.py sum
```

## Ongoing Development

This code is being developed on an on-going basis at the author's
[Github site](https://github.com/tkralphs/JoCTemplate).

## Support

For support in using this software, submit an
[issue](https://github.com/tkralphs/JoCTemplate/issues/new).
