
# A Modified Sequential Probability Ratio Test (MSPRT)

Article: [Pramanik, Johnson and Bhattacharya (2020+)](https://arxiv.org/abs/1811.08478)


## Description

Given the maximum available sample size $(N)$ for an experiment, and the target levels of Type I and II error probabilities, this package designs a modified SPRT (MSPRT). For any designed MSPRT the package can also obtain its operating characteristics and implement the test for a given sequentially observed data. The MSPRT is defined in a manner very similar to Wald's initial proposal. The proposed test has shown evidence of reducing the average sample size required to perform statistical hypothesis tests at specified levels of significance and power. Currently, the package implements one-sample proportion tests, one and two-sample $z$ tests, and one and two-sample $t$ tests. A brief user guidance for this package is provided below. One can also refer to the supplemental information for the same.
             

## Key functions

The followings are the key functions in the package:

  * [**Type2.fixed_design**](#type2.fixed_design)
  * [**fixed_design.alt**](#fixed_design.alt)
  * [**UMPBT.alt**](#umpbt.alt)
  * [**design.MSPRT**](#design.msprt)
  * [**OCandASN.MSPRT**](#ocandasn.msprt)
  * [**implement.MSPRT**](#implement.msprt)
  * [**effectiveN.oneProp**](#effectiven.oneprop)
  * [**Nstar**](#nstar)
             

## User's guide to the key functions


### **Type2.fixed_design**

<font size="3"> [[List of all functions]](#key-functions) </font>

Obtains the Type II error probability of fixed-design tests for testing the point null hypothesis $H_0 : \theta = \theta_0$.

* **Input:**

  * **theta:** Numeric. Effect size where the Type II error probability is desired.
  
  * **test.type:** Character. Type of test. Currently, the package only allows 
    * 'oneProp' for one-sample proportion tests
    * 'oneZ' for one-sample $z$ tests
    * 'oneT' for one-sample $t$ tests
    * 'twoZ' for two-sample $z$ tests
    * 'twoT' for two-sample $t$ tests.
    
  * **side:** Character. Direction of the composite alternative hypothesis. 'right' for $H_1 : \theta > \theta_0$ (**default**), and 'left' for $H_1 : \theta < \theta_0$.
  
  * **theta0:** Numeric. Hypothesized value of effect size ($\theta_0$) under $H_0$. **Default:** 0.5 in one-sample proportion tests, and 0 for others.
  
  * **N:** Positive integer. Sample size in one-sample tests.
  
  * **N1:** Positive integer. Sample size from Group-1 in two-sample tests.
  
  * **N2:** Positive integer. Sample size from Group-2 in two-sample tests.
  
  * **Type1:** Numeric in $[0,1]$. Prespecified Type I error probability. **Default:** 0.005.
  
  * **sigma:** Positive numeric. Known standard deviation in one-sample $z$ tests. **Default:** 1.
  
  * **sigma1:** Positive numeric. Known standard deviation for Group-1 in two-sample $z$ tests. **Default:** 1.
  
  * **sigma2:** Positive numeric. Known standard deviation for Group-2 in two-sample $z$ tests. **Default:** 1.

* **Output:** Numeric in $[0,1]$. The Type II error probability of the fixed-design test at the specified effect size value 'theta'.


### **fixed_design.alt** 

<font size="3"> [[List of all functions]](#key-functions) </font>

Given a sample size and prespecified Type I & II error probabilities, this function obtains the fixed-design alternative ($\theta_a$) for testing the point null hypothesis $H_0 : \theta = \theta_0$. 

**Note:** 'fixed-design alternative' is the effect size under $H_1$ where the fixed design test with the prespecified Type I error probability and the given sample size has the prespecified Type II error probability.

* **Input:**

  * **test.type:** Character. Type of test. Currently, the package only allows 
    * 'oneProp' for one-sample proportion tests
    * 'oneZ' for one-sample $z$ tests
    * 'oneT' for one-sample $t$ tests
    * 'twoZ' for two-sample $z$ tests
    * 'twoT' for two-sample $t$ tests.
    
  * **side:** Character. Direction of the composite alternative hypothesis. 'right' for $H_1 : \theta > \theta_0$ (**default**), and 'left' for $H_1 : \theta < \theta_0$.
  
  * **theta0:** Numeric. Hypothesized value of effect size ($\theta_0$) under $H_0$. **Default:** 0.5 in one-sample proportion tests, and 0 for others.
  
  * **N:** Positive integer. Sample size in one-sample tests.
  
  * **N1:** Positive integer. Sample size from Group-1 in two-sample tests.
  
  * **N2:** Positive integer. Sample size from Group-2 in two-sample tests.
  
  * **Type1:** Numeric within $[0,1]$. Prespecified Type I error probability. **Default:** 0.005.
  
  * **Type2:** Numeric in $[0,1]$. Prespecified Type II error probability. **Default:** 0.2.
  
    * **sigma:** Positive numeric. Known standard deviation in one-sample $z$ tests. **Default:** 1.
    
  * **sigma1:** Positive numeric. Known standard deviation for Group-1 in two-sample $z$ tests. **Default:** 1.
  
  * **sigma2:** Positive numeric. Known standard deviation for Group-2 in two-sample $z$ tests. **Default:** 1.

* **Output:** Numeric. The fixed-design alternative effect size ($\theta_a$).


### **UMPBT.alt** 

<font size="3"> [[List of all functions]](#key-functions) </font>

Given a sample size and prespecified Type I & II error probabilities, this function obtains the objective alternative in the Uniformly Most Powerful Bayesian Test (UMPBT) [[Johnson (2013a)](https://projecteuclid.org/euclid.aos/1378386237#info), [Johnson (2013b)](https://www.stat.tamu.edu/~vjohnson/files/PNAS2013.pdf)].

* **Input:**

  * **test.type:** Character. Type of test. Currently, the package only allows 
    * 'oneProp' for one-sample proportion tests
    * 'oneZ' for one-sample $z$ tests
    * 'oneT' for one-sample $t$ tests
    * 'twoZ' for two-sample $z$ tests
    * 'twoT' for two-sample $t$ tests.
    
  * **side:** Character. Direction of the composite alternative hypothesis. 'right' for $H_1 : \theta > \theta_0$ (**default**), and 'left' for $H_1 : \theta < \theta_0$.
  
  * **theta0:** Numeric. Hypothesized value of effect size ($\theta_0$) under $H_0$. **Default:** 0.5 in one-sample proportion tests, and 0 for others.
  
  * **N:** Positive integer. Sample size in one-sample tests.
  
  * **N1:** Positive integer. Sample size from Group-1 in two-sample tests.
  
  * **N2:** Positive integer. Sample size from Group-2 in two-sample tests.
  
  * **Type1:** Numeric within $[0,1]$. Prespecified Type I error probability. **Default:** 0.005.
  
  * **sigma:** Positive numeric. Known standard deviation in one-sample $z$ tests. **Default:** 1.
  
  * **sigma1:** Positive numeric. Known standard deviation for Group-1 in two-sample $z$ tests. **Default:** 1.
  
  * **sigma2:** Positive numeric. Known standard deviation for Group-2 in two-sample $z$ tests. **Default:** 1.
  
  * **obs:** Numeric vector. The vector of observations based on which the UMPBT alternative in one-sample $t$ test is determined. *Either 'obs' or 'sd.obs' is required*.
  
  * **sd.obs:** Positive numeric. The standard deviation (with divisor $n-1$) of observations based on which the UMPBT alternative in one-sample $t$ test is determined. *Either 'obs' or 'sd.obs' is required.*
  
  * **obs1:** Numeric vector. The vector of observations from Group-1 based on which the UMPBT alternative in two-sample $t$ test is determined. *Either both 'obs1' and 'obs2', or 'pooled.sd' is required.*
  
  * **obs2:** Numeric vector. The vector of observations from Group-2 based on which the UMPBT alternative in two-sample $t$ test is determined. *Either both 'obs1' and 'obs2', or 'pooled.sd' is required.*
  
  * **pooled.sd:** Positive numeric. The pooled standard deviation of observations from Group-1 and 2 based on which the UMPBT alternative in two-sample $t$ test is determined. *Either both 'obs1' and 'obs2', or 'pooled.sd' is required.*

* **Output:** List with two named components 'theta' and 'mix.prob' in one-sample proportion test. In this case, the UMPBT alternative is a mixture distribution of two points. 'theta' contains the two points (effect sizes) and 'mix.prob' contains their respective mixing probabilities.

  Numeric in case of all the other tests. It is the UMPBT alternative effect size.


### **design.MSPRT** 

<font size="3"> [[List of all functions]](#key-functions) </font>

Given the maximum available sample size and prespecified Type I & II error probabilities, this function designs/obtains the corresponding MSPRT [[Pramanik, Johnson and Bhattacharya (2020+)](https://arxiv.org/abs/1811.08478)].

* **Input:**

  * **test.type:** Character. Type of test. Currently, the package only allows 
    * 'oneProp' for one-sample proportion tests
    * 'oneZ' for one-sample $z$ tests
    * 'oneT' for one-sample $t$ tests
    * 'twoZ' for two-sample $z$ tests
    * 'twoT' for two-sample $t$ tests.
    
  * **side:** Character. Direction of the composite alternative hypothesis. 
    * 'right' to test against the right-sided alternative $H_1 : \theta > \theta_0$ (**default**).
    * 'left' to test against the left-sided alternative $H_1 : \theta < \theta_0$.
    * 'both' to test against the two-sided alternative $H_1 : \theta \neq \theta_0$.
  
  * **theta0:** Numeric. Hypothesized value of effect size ($\theta_0$) under $H_0$. **Default:** 0.5 in one-sample proportion tests, and 0 for others.
  
  * **theta1:** Logical, numeric or list (two components with names 'right' and 'left').
  
    * If 'FALSE', no comparison is done under the alternative hypothesis.
    * If 'TRUE' (**Default**), comparison is done at the fixed-design alternative effect size ($\theta_a$).
    * If numeric, this can only be in case of one-sided tests (that is, side = 'right' or 'left'). The comparison is done at the specified numeric value of the alternative effect size.
    * If list, this can only be in case of two-sided tests (that is, side = 'both'). The list has to be of the form list('right' = $\theta_1$, 'left' = $\theta_2$). Then the comparison is done at alternative effect sizes $\theta_1$ and $\theta_2$. 
    
    <u>Note:</u> In case of two-sided tests at a given level of significance, there are two effect sizes under $H_1$ (one on the right of $H_0$ and one on the left) that corresponds to the same Type II error probability (or power). This list provides users with the ability where he/she can replace $\theta_1$ and $\theta_2$ by any effect sizes from each side in the form of a list as mentioned above, and can get the designed MSPRT together with its operating characteristics at those effect sizes.
    
  * **Type1.target:** Numeric within $[0,1]$. Prespecified level of Type I error probability. **Default:** 0.005. The MSPRT exactly maintains its Type I error probability at this value.
  
  * **Type2.target:** Numeric within $[0,1]$. Prespecified level of Type 2 error probability. **Default:** 0.2. The MSPRT approximately maintains its Type II error probability at this value at the corresponding fixed-design alternative ($\theta_a$).
  
  * **N.max:** Positive integer. Maximum available sample size in one-sample tests.
  
  * **N1.max:** Positive integer. Maximum available sample size from Group-1 in two-sample tests.
  
  * **N2.max:** Positive integer. Maximum available sample size from Group-2 in two-sample tests.
  
  * **sigma:** Positive numeric. Known standard deviation in one-sample $z$ tests. **Default:** 1.
  
  * **sigma1:** Positive numeric. Known standard deviation for Group-1 in two-sample $z$ tests. **Default:** 1.
  
  * **sigma2:** Positive numeric. Known standard deviation for Group-2 in two-sample $z$ tests. **Default:** 1.
  
  * **batch.size:** Integer vector. A vector denoting the number of observations that are planned to be observed at each sequential step in one-sample tests. 
  
    * **Default:** 
      * <u>Proportion and $z$ tests:</u> rep(1, *N.max*).
      * <u>$t$ tests:</u> c(2, rep(1, *N.max*-1)).
      
    Default values mean the sequential analysis is performed after observing each observation. This corresponds to a sequential MSPRT. If any batch size is more than 1 (or more than 2 in the 1st step for $t$ test) it corresponds to a group sequential MSPRT.
    <u> Note:</u> First batch size for $t$ tests needs to be at least 2. The length of batch.size equals to the maximum number of planned sequential analyses.
    
  * **batch1.size:** Integer vector. A vector denoting the number of observations that are planned to be observed from Group-1 at each sequential step in two-sample tests. 
    * **Default:** 
      * <u> $z$ tests:</u> rep(1, *N1.max*).
      * <u> $t$ tests:</u> c(2, rep(1, *N1.max*-1)).
    
    Default values mean the sequential analysis is performed after observing each observation from Group-1.
  * **batch2.size:** Integer vector. A vector denoting the number of observations that are planned to be observed from Group-2 at each sequential step in two-sample tests. 
    * **Default:** 
      * <u> $z$ tests:</u> rep(1, *N2.max*).
      * <u> $t$ tests:</u> c(2, rep(1, *N2.max*-1)).
    
    Default values mean the sequential analysis is performed after observing each observation from Group-2. 
    
    <u> Note:</u> The default values of batch1.size and batch2.size correspond to a sequential MSPRT. If any of batch1.size or batch2.size is more than 1 (or more than 2 in the 1st step for $t$ test) it corresponds to a group sequential MSPRT. First batch size for $t$ tests needs to be at least 2. The length of batch1.size and batch2.size should be equal, and it is the maximum number of planned sequential analyses.
    
  * **nReplicate:** Positive integer. Total number of replications to be used in Monte Carlo simulation for calculating the termination threshold and the operating characteristics of the MSPRT. **Default:** $10^6$.
    
  * **verbose:** Logical. If TRUE (**default**), returns messages of the current proceedings. Otherwise it doesn't.
  
  * **seed:** Integer. Random number generating seed. **Default:** $1$.
    

* **Output:** List. The list has the following named components in case of one-sided one-sample tests:

  * **TypeI.attained:** Numeric in $[0,1]$. Type I error probability attained by the designed MSPRT.
  
  * **Type2.attained:** Numeric in $[0,1]$. Type II error probability attained by the designed MSPRT at the specified alternative effect size *theta1*. Returned only if *theta1* is TRUE or numeric.
  
  * **N:** List.
    * If *theta1* = FALSE, the list has one component named *H0*. It stores an integer vector of length *nReplicate*. This is the vector of sample size required by the MSPRT for each of *nReplicate* Monte Carlo simulations under $H_0$.
    * If *theta1* is TRUE or numeric, the list has two components named *H0* and *H1*. Each of these stores an integer vector of length *nReplicate*. The stored vector under *H0* is the same as in *theta1* = FALSE. The *H1* component stored the vector of sample size required by the MSPRT for each of *nReplicate* Monte Carlo simulations under the specified alternative effect size.
      
  * **EN:** Numeric vector.
      * If *theta1* = FALSE, the vector is of length 1. It is the number of samples required on average by the MSPRT under $H_0$.
      * If theta1 is TRUE or numeric, the vector is of length 2. They are the number of samples required on average by the MSPRT under $H_0$ (first component) and the specified alternative effect size (second component), respectively.
      
  * **UMPBT** or **theta.UMPBT:** The UMPBT alternative. *UMPBT* in case of one-sample proportion test and *theta.UMPBT* in case of all the other tests. Their types are the same as their output from *UMPBT.alt* function.
  
    <u> Note:</u> Not returned in $t$ tests as it depends on the data.
    
  * **theta1:** Numeric. Returned only if *theta1* is anything but FALSE. Stores the effect size under $H_1$ where the operating characteristic of the MSPRT is obtained.
  
  * **Type2.fixed.design:** Numeric in $[0,1]$. Type II error probability attained by the fixed design test with sample size *N.max* and Type I error probability *Type1.target* at the alternative effect size *theta1*.
  
  * **RejectH0.threshold:** Positive numeric. Threshold for rejecting $H_0$ in the MSPRT.
  
  * **RejectH1.threshold:** Positive numeric. Threshold for accepting $H_1$ in the MSPRT.
  
  * **termination.threshold:** Positive numeric. Termination threshold of the MSPRT.
  
  In case of one-sided two-sample tests the above components are returned with following modifications:

    * **N:**  List.
      * If *theta1* = FALSE the list has one component named *H0*.
      * If *theta1* is TRUE or numeric, the list has two components named *H0* and *H1*.
    
      Each of the named components *H0* and *H1* contains a list with two components named *Group1* and *Group2*. Each of these contains the same vector corresponding to Group-1 and Group-2. In each of these, it contains the sample size required by the MSPRT in each of *nReplicate* Monte Carlo simulations under the respective effect size for the respective group.
  
  * **EN:** List.
    * If *theta1* = FALSE the list has one component named *H0*.
    * If *theta1* is TRUE or numeric, the list has two components named *H0* and *H1*.
    
    Each of the named components *H0* or *H1* contains a list with two components named *Group1* and *Group2*. In each of these, it contains the sample size required on average by the MSPRT under the respective effect size for the respective group.
  
  In case of two-sided tests the above components are returned with following modifications:

    * **Type2.attained:** Numeric vector of length 2 with both elements in $[0,1]$. The first and second component is the Type II error probability of the MSPRT at the specified alternative effect sizes *theta1$\$$right* and *theta1$\$$left*, respectively.
  
    * **N:** This is the same as in one-sided tests if *theta1* is FALSE. If *theta1* is TRUE or a two-component list with names *right* and *left*, this is a list with three components with names *H0*, *right* and *left* instead of a two-component list with names *H0* and *H1*. Quantities stored under these components are the same as in one-sided tests except the quantities under *right* and *left* are the same performance of the designed MSPRT at the specified alternative effect sizes *theta1$\$$right* and *theta1$\$$left*, respectively.
  
    * **EN:** Numeric vector. The same as in one-sided tests if *theta1* is FALSE. If *theta1* is TRUE or a two-component list with names *right* and *left*, this is a numeric vector of length 3, where the first, second and third components are the average required sample size under $H_0$, and at the specified alternative effect sizes *theta1$\$$right* and *theta1$\$$left*, respectively.
  
  Additionally, the output list also contains the provided arguments of design.MSPRT, and 
  
    * **nAnalyses:** Positive integer. This is the maximum number of sequential analyses that is planned. This equals to the *length(batch.size)* in one-sample tests, and to the *length(batch1.size)* and *length(batch2.size)* in two-sample tests.


### **OCandASN.MSPRT** 

<font size="3"> [[List of all functions]](#key-functions) </font>

This function obtains the operating characteristics, that is the probability of accepting $H_0$ and the sample size required on average for reaching a decision, for a designed MSPRT at the specified effect size(s).

* **Input:**

  * **theta:** Numeric vector. Vector of effect size(s) where the operating characteristics of the MSPRT is desired.
  
  * **design.MSPRT.object:** List. The output returned from [design.MSPRT](#design.msprt) corresponding to the MSPRT for which the operating characteristics are desired.
    
  * **nReplicate:** Positive integer. Total number of replications to be used in Monte Carlo simulation for calculating the operating characteristics of the MSPRT. **Default:** $10^6$.
    
  * **verbose:** Logical. If TRUE (**default**), returns messages of the current proceedings. Otherwise it doesn't.
  
  * **nCore:** Positive integer. Total number of cores available for computation. Can be anything $\geq 1$. **Default:** *detectCores( ) - 1*. That is, 1 less than the total number of available cores.
  
  * **seed:** Positive integer. Random number generating seed. **Default:** $1$.
  
  The rest of the arguments are similar to that of [design.MSPRT](#design.msprt). If *design.MSPRT.object* is provided, no other arguments are required (Easier option). Otherwise, they are required just like in [design.MSPRT](#design.msprt), except
  
  * **termination.threshold:** Positive numeric. Termination threshold of the MSPRT.


* **Output:** Data frame.
  * <u>One-sample tests:</u> The data frame has 3 columns named *theta*, *acceptH0.prob* and *EN*, and the number of rows equals to the number of effect sizes (length of *theta*) where the operating characteristics are evaluated. Each row corresponds to a particular value of theta (effect size). The columns respectively contain the value of a particular theta (effect size), and the probability of accepting the $H_0$ and the average sample size required by the MSPRT for reaching a decision thereat.
  
  * <u>Two-sample tests:</u> The data frame has 4 columns named *theta*, *acceptH0.prob*, *EN1* and *EN2*, and the number of rows equals to the number of effect sizes (length of *theta*) where the operating characteristics are evaluated. Each row corresponds to a particular value of theta (effect size). The columns respectively contain the value of a particular theta (effect size), and the probability of accepting the $H_0$ at that effect size, and the average sample size from Group-1 & 2 that is required by the MSPRT for reaching a decision thereat.


### **implement.MSPRT** 

<font size="3"> [[List of all functions]](#key-functions) </font>

This function implements the MSPRT for a sequentially observed data.

* **Input:**

  * **obs:** Numeric vector. The vector of data in the order they are sequentially observed for one-sample tests.
  <u>Note:</u> Its length
  
  * **obs1:** Numeric vector. The vector of data in the order they are sequentially observed from Group-1 for two-sample tests.
  
  * **obs2:** Numeric vector. The vector of data in the order they are sequentially observed from Group-2 for two-sample tests.
    
  * **verbose:** Logical. If TRUE (**default**), returns messages of the current proceedings. Otherwise it doesn't.
    
  * **plot.it:** Logical. If TRUE (**default**), returns a sequential comparison plot. Otherwise it doesn't.
    

* **Output:** List. The list has the following named components in case of one-sided one-sample tests:

  * **n:**  Positive integer. Number of samples required to reach the decision.
  
  * **decision:**  Character. The decision reached. The possibilities are 'accept', 'reject' and 'continue'. They respectively correspond to accepting $H_0$, rejecting $H_0$ and continue sampling.
  
  * **RejectH0.threshold:** Positive numeric. Threshold for rejecting $H_0$ in the MSPRT.
  
  * **RejectH1.threshold:** Positive numeric. Threshold for accepting $H_1$ in the MSPRT.
  
  * **LR:** Numeric vector. Vector of weighted likelihood ratios (proportion tests) or likelihood ratios ($z$ tests) or Bayes factor ($t$ tests) that are computed at each step of sequential analysis until either a decision is reached or the maximum available number of samples (*N.max* in one-sample tests, or *N1.max* and *N2.max* in two-sample tests) has been used.
  
  * **UMPBT alternative** This stores the UMPBT alternative(s) as
    * **UMPBT** for proportion tests. Of the same type as it is returned by [UMPBT.alt](#umpbt.alt) in these tests.
    * **theta.UMPBT** for $z$ and $t$ tests. This is a numeric in case of $z$ tests and a numeric vector in case of $t$ tests. For $t$ tests the UMPBT alternative depends on the data. So the numeric vector returned in this case contains the UMPBT alternative computed at step of sequential analysis and is based on all data observed until that step.
    
  In case of two-sample tests, the *n* output above is replaced by *n1* and *n2*. They are positive integers and refer to the number of samples from Group-1 and 2 required to reach the decision.
  
  In case of two-sided tests at level of significance $\alpha$, the MSPRT carries out a right and a left sided test simultaneously at level of significance $\alpha/2$. In this case the outputs are same as above with following changes in components in the returned list:
  
    * **LR** List. It has two components named *'right'* and *'left'* corresponding to the right and left sided tests of size $\alpha/2$. Each of these components stores the vector of weighted likelihood ratios (proportion tests) or likelihood ratios ($z$ tests) or Bayes factor ($t$ tests) that are computed at each step of sequential analysis until either a decision is reached or the maximum available number of samples (*N.max* in one-sample tests, or *N1.max* and *N2.max* in two-sample tests) has been used for that sided test.
    
    * *UMPBT* or *theta.UMPBT* is a list with two components named *'right'* and *'left'* corresponding to the right and left sided tests of size $\alpha/2$. Each of these contains the UMPBT alternative (of the same type as the output from [UMPBT.alt](#umpbt.alt)) for the test with respective sides.


### **effectiveN.oneProp** 

<font size="3"> [[List of all functions]](#key-functions) </font>

Given a maximum sample size that is planned to use, this function obtains the maximum sample size ($N$) that is suggested to use in designing the MSPRT for one-sample proportion tests.

* **Input:**
  
  * **N:** Positive integer. Maximum sample that is intended to use.
    
  * **side:** Character. Direction of the composite alternative hypothesis. 'right' for $H_1 : \theta > \theta_0$ (**default**), and 'left' for $H_1 : \theta < \theta_0$.
  
  * **Type1:** Numeric in $[0,1]$. Prespecified Type I error probability. **Default:** 0.005.
  
  * **theta0:** Numeric. Hypothesized value of effect size ($\theta_0$) under $H_0$. **Default:** 0.5.
    
  * **plot.it:** Logical. If TRUE (**default**), returns a plot. Otherwise it doesn't.

* **Output:** Positive integer. This is suggested to use in [design.MSPRT](#design.msprt) as the maximum availeble sample size($N$) to design the MSPRT for one-sample proportion tests.


### **Nstar** 

<font size="3"> [[List of all functions]](#key-functions) </font>

Given the sample size that is available at a lower level of significance, this function calculates the sample size that is required for achieving a higher level of significance so that a desired level of Type II error probability is maintained at a desired effect size.

* **Input:**

  * **test.type:** Character. Type of test. Currently, the package only allows 
    * 'oneProp' for one-sample proportion tests
    * 'oneZ' for one-sample $z$ tests
    * 'oneT' for one-sample $t$ tests
    * 'twoZ' for two-sample $z$ tests
    * 'twoT' for two-sample $t$ tests.
  
  * **N:** Positive integer. Sample size available at the lower level of significance in one-sample tests.
  
  * **N1:** Positive integer. Sample size available from Group-1 at the lower level of significance in two-sample tests.
  
  * **N2:** Positive integer. Sample size available from Group-2 at the lower level of significance in two-sample tests.
  
  * **N.increment:** Positive integer. Increment in sample size allowed while searching for the sample size that is required for achieving the higher level of significance.
  
  * **N1.increment:** Positive integer. Increment in sample size from Group-1 allowed while searching for the sample size that is required for achieving the higher level of significance.
  
  * **N1.increment:** Positive integer. Increment in sample size from Group-2 allowed while searching for the sample size that is required for achieving the higher level of significance.
  
  * **lower.signif:** Numeric within $[0,1]$. Lower level of significance. **Default:** 0.05.
  
  * **higher.signif:** Numeric within $[0,1]$. Higher level of significance. **Default:** 0.005.
  
  * **theta0:** Numeric. Hypothesized value of effect size ($\theta_0$) under $H_0$. **Default:** 0.5 in one-sample proportion tests, and 0 for others.
  
  * **side:** Character. Direction of the composite alternative hypothesis. 'right' for $H_1 : \theta > \theta_0$ (**default**), and 'left' for $H_1 : \theta < \theta_0$.
  
  * **Type2.target:** Numeric in $[0,1]$. Desired Type II error probability at effect size *theta*. **Default:** 0.2.
  
  * **theta:** Numeric. Effect size value where *Type2.target* Type II error probability is desired at both levels of significance. **Default:** Fixed-design alternative ($\theta_a$) at the lower level of significance; that is, the effect size where the fixed design test with $N$ samples and level of significance *lower.signif* has the Type II error probability *Type2.target*.
  
  * **sigma:** Positive numeric. Known standard deviation in one-sample $z$ tests. **Default:** 1.
    
  * **sigma1:** Positive numeric. Known standard deviation for Group-1 in two-sample $z$ tests. **Default:** 1.
  
  * **sigma2:** Positive numeric. Known standard deviation for Group-2 in two-sample $z$ tests. **Default:** 1.
    
  * **plot.it:** Logical. If TRUE (**default**), returns a plot. Otherwise it doesn't.

* **Output:** 

  * <u>One-sample tests:</u> Numeric. The required sample size.
  
  * <u>Two-sample tests:</u> Numeric vector of length 2. The first and second components store the sample sizes required respectively from Group 1 and 2 for achieving the higher level of significance.




