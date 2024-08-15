# Changepoint Detection

## Overview

This repository contains implementations of changepoint detection algorithms using two different approaches:

1. **Fixed Number of Segments (k)**: Detecting changepoints by specifying the desired number of segments.
2. **Penalty Parameter (λ)**: Detecting changepoints by controlling the penalty parameter, which influences the number of changepoints detected.

These methods are useful for analyzing time series data to identify points where the statistical properties of the data change significantly.

## Table of Contents

- [Introduction](#introduction)
- [License](#license)

## Introduction

Changepoint detection is an important technique in time series analysis, helping to identify moments where the data behavior changes, such as shifts in mean, variance, or other statistical properties. This repository provides two methods for detecting changepoints:

1. **Fixed Number of Segments (k)**: This method allows the user to specify the number of segments they believe are present in the data, and the algorithm will detect changepoints accordingly.
   
2. **Penalty Parameter (λ)**: This method allows the user to control the trade-off between the number of changepoints detected and the fit of the model. A higher penalty will result in fewer changepoints, while a lower penalty will detect more changepoints.

## License
This project is licensed under the MIT License. See the LICENSE file for more details.
