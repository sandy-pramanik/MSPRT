
source('MSPRT.R')

#################### one-sample proportion test ####################
N.max = 20

#### right sided ####
# design
design.oneprop.right = design.MSPRT(test.type = 'oneProp', side = 'right',
                                    N.max = N.max)

# OC and ASN
OC.oneprop.right = OCandASN.MSPRT(theta = seq(design.oneprop.right$theta0, 1,
                                              length.out = 5),
                                  design.MSPRT.object = design.oneprop.right)

# implementation
set.seed(1)
theta.gen = 0.5  # change effect size to experiment
y = rbinom(N.max, 1, theta.gen)
implement.oneprop.right = implement.MSPRT(obs = y, 
                                          design.MSPRT.object = design.oneprop.right)

#### left sided ####
# design
design.oneprop.left = design.MSPRT(test.type = 'oneProp', side = 'left',
                                    N.max = N.max)

# OC and ASN
OC.oneprop.left = OCandASN.MSPRT(theta = seq(0, design.oneprop.right$theta0,
                                             length.out = 5),
                                  design.MSPRT.object = design.oneprop.left)

# implementation
set.seed(1)
theta.gen = 0.5  # change effect size to experiment
y = rbinom(N.max, 1, theta.gen)
implement.oneprop.left = implement.MSPRT(obs = y, 
                                          design.MSPRT.object = design.oneprop.left)

#### both sided ####
# design
design.oneprop.both = design.MSPRT(test.type = 'oneProp', side = 'both',
                                   N.max = N.max)

# OC and ASN
OC.oneprop.both = OCandASN.MSPRT(theta = seq(0, 1, length.out = 5),
                                 design.MSPRT.object = design.oneprop.both)

# implementation
set.seed(1)
theta.gen = 0.5  # change effect size to experiment
y = rbinom(N.max, 1, theta.gen)
implement.oneprop.both = implement.MSPRT(obs = y, 
                                         design.MSPRT.object = design.oneprop.both)


#################### one-sample z test ####################
N.max = 20

#### right sided ####
# design
design.onez.right = design.MSPRT(test.type = 'oneZ', side = 'right',
                                    N.max = N.max)

# OC and ASN
OC.onez.right = OCandASN.MSPRT(theta = seq(design.onez.right$theta0,
                                           design.onez.right$theta0 + 3*design.onez.right$sigma,
                                           length.out = 5),
                               design.MSPRT.object = design.onez.right)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y = rnorm(N.max, theta.gen, design.onez.right$sigma)
implement.onez.right = implement.MSPRT(obs = y, 
                                       design.MSPRT.object = design.onez.right)

#### left sided ####
# design
design.onez.left = design.MSPRT(test.type = 'oneZ', side = 'left',
                                N.max = N.max)

# OC and ASN
OC.onez.left = OCandASN.MSPRT(theta = seq(design.onez.left$theta0 - 3*design.onez.left$sigma,
                                          design.onez.left$theta0,
                                          length.out = 5),
                              design.MSPRT.object = design.onez.left)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y = rnorm(N.max, theta.gen, design.onez.left$sigma)
implement.onez.left = implement.MSPRT(obs = y, 
                                      design.MSPRT.object = design.onez.left)

#### both sided ####
# design
design.onez.both = design.MSPRT(test.type = 'oneZ', side = 'both',
                                N.max = N.max)

# OC and ASN
OC.onez.both = OCandASN.MSPRT(theta = seq(design.onez.both$theta0 - 3*design.onez.both$sigma,
                                          design.onez.both$theta0 + 3*design.onez.both$sigma,
                                          length.out = 5),
                              design.MSPRT.object = design.onez.both)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y = rnorm(N.max, theta.gen, design.onez.both$sigma)
implement.onez.both = implement.MSPRT(obs = y, 
                                      design.MSPRT.object = design.onez.both)


#################### one-sample t test ####################
N.max = 20

#### right sided ####
# design
design.onet.right = design.MSPRT(test.type = 'oneT', side = 'right',
                                 N.max = N.max)

# OC and ASN
OC.onet.right = OCandASN.MSPRT(theta = seq(design.onet.right$theta0, 1,
                                           length.out = 5),
                               design.MSPRT.object = design.onet.right)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y = rnorm(N.max, theta.gen, 1)
implement.onet.right = implement.MSPRT(obs = y, 
                                       design.MSPRT.object = design.onet.right)

#### left sided ####
# design
design.onet.left = design.MSPRT(test.type = 'oneT', side = 'left',
                                N.max = N.max)

# OC and ASN
OC.onet.left = OCandASN.MSPRT(theta = seq(-1, design.onet.left$theta0,
                                          length.out = 5),
                              design.MSPRT.object = design.onet.left)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y = rnorm(N.max, theta.gen, 1)
implement.onet.left = implement.MSPRT(obs = y, 
                                      design.MSPRT.object = design.onet.left)

#### both sided ####
# design
design.onet.both = design.MSPRT(test.type = 'oneT', side = 'both',
                                N.max = N.max)

# OC and ASN
OC.onet.both = OCandASN.MSPRT(theta = seq(-1, 1, length.out = 5),
                              design.MSPRT.object = design.onet.both)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y = rnorm(N.max, theta.gen, 1)
implement.onet.both = implement.MSPRT(obs = y, 
                                      design.MSPRT.object = design.onet.both)


#################### two-sample z test ####################
N1.max = N2.max = 20

#### right sided ####
# design
design.twoz.right = design.MSPRT(test.type = 'twoZ', side = 'right',
                                 N1.max = N1.max, N2.max = N2.max)

# OC and ASN
OC.twoz.right = OCandASN.MSPRT(theta = seq(design.twoz.right$theta0, 
                                           design.twoz.right$theta0 + 2,
                                           length.out = 5),
                               design.MSPRT.object = design.twoz.right)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y1 = rnorm(N1.max, theta.gen/2, design.twoz.right$sigma1)
y2 = rnorm(N2.max, -theta.gen/2, design.twoz.right$sigma2)
implement.twoz.right = implement.MSPRT(obs1 = y1, obs2 = y2,
                                       design.MSPRT.object = design.twoz.right)

#### left sided ####
# design
design.twoz.left = design.MSPRT(test.type = 'twoZ', side = 'left',
                                N1.max = N1.max, N2.max = N2.max)

# OC and ASN
OC.twoz.left = OCandASN.MSPRT(theta = seq(design.twoz.left$theta0 - 2,
                                          design.twoz.left$theta0,
                                          length.out = 5),
                              design.MSPRT.object = design.twoz.left)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y1 = rnorm(N1.max, theta.gen/2, design.twoz.left$sigma1)
y2 = rnorm(N2.max, -theta.gen/2, design.twoz.left$sigma2)
implement.twoz.left = implement.MSPRT(obs1 = y1, obs2 = y2,
                                      design.MSPRT.object = design.twoz.left)

#### both sided ####
# design
design.twoz.both = design.MSPRT(test.type = 'twoZ', side = 'both',
                                N1.max = N1.max, N2.max = N2.max)

# OC and ASN
OC.twoz.both = OCandASN.MSPRT(theta = seq(design.twoz.both$theta0 - 2,
                                          design.twoz.both$theta0 + 2,
                                          length.out = 5),
                              design.MSPRT.object = design.twoz.both)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y1 = rnorm(N1.max, theta.gen/2, design.twoz.both$sigma1)
y2 = rnorm(N2.max, -theta.gen/2, design.twoz.both$sigma2)
implement.twoz.both = implement.MSPRT(obs1 = y1, obs2 = y2,
                                      design.MSPRT.object = design.twoz.both)


#################### two-sample t test ####################
N1.max = N2.max = 20

#### right sided ####
# design
design.twot.right = design.MSPRT(test.type = 'twoT', side = 'right',
                                 N1.max = N1.max, N2.max = N2.max)

# OC and ASN
OC.twot.right = OCandASN.MSPRT(theta = seq(design.twot.right$theta0, 
                                           design.twot.right$theta0 + 2,
                                           length.out = 5),
                               design.MSPRT.object = design.twot.right)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y1 = rnorm(N1.max, theta.gen/2, 1)
y2 = rnorm(N2.max, -theta.gen/2, 1)
implement.twot.right = implement.MSPRT(obs1 = y1, obs2 = y2,
                                       design.MSPRT.object = design.twot.right)

#### left sided ####
# design
design.twot.left = design.MSPRT(test.type = 'twoT', side = 'left',
                                N1.max = N1.max, N2.max = N2.max)

# OC and ASN
OC.twot.left = OCandASN.MSPRT(theta = seq(design.twot.left$theta0 - 2,
                                          design.twot.left$theta0,
                                          length.out = 5),
                              design.MSPRT.object = design.twot.left)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y1 = rnorm(N1.max, theta.gen/2, 1)
y2 = rnorm(N2.max, -theta.gen/2, 1)
implement.twot.left = implement.MSPRT(obs1 = y1, obs2 = y2,
                                      design.MSPRT.object = design.twot.left)

#### both sided ####
# design
design.twot.both = design.MSPRT(test.type = 'twoT', side = 'both',
                                N1.max = N1.max, N2.max = N2.max)

# OC and ASN
OC.twot.both = OCandASN.MSPRT(theta = seq(design.twot.both$theta0 - 2,
                                          design.twot.both$theta0 + 2,
                                          length.out = 5),
                              design.MSPRT.object = design.twot.both)

# implementation
set.seed(1)
theta.gen = 0  # change effect size to experiment
y1 = rnorm(N1.max, theta.gen/2, 1)
y2 = rnorm(N2.max, -theta.gen/2, 1)
implement.twot.both = implement.MSPRT(obs1 = y1, obs2 = y2,
                                      design.MSPRT.object = design.twot.both)


