# Code adapted from Wikipedia

####################################################
############  START MODIFICATIONS HERE  ############
####################################################
'''
states = ('Fair', 'Biased')
 
#observations = ('tails', 'heads', 'heads', 'heads', 'heads')
observations = ('tails', 'tails', 'heads', 'heads', 'heads')

start_probability = {'Fair': 0.50, 'Biased': 0.50}
 
transition_probability = {
   'Fair' : {'Fair': 0.90, 'Biased': 0.10},
   'Biased' : {'Fair': 0.10, 'Biased': 0.90}
   }
 
emission_probability = {
   'Fair' : {'heads': 0.50, 'tails': 0.50},
   'Biased' : {'heads': 0.75, 'tails': 0.25}
   }

states = ('exon', 'donor1', 'donor2', 'intron', 'acceptor1', 'acceptor2')
 
observations = ('A','T','G','G','C','C','C','G','A','A','C','C','A','A','G','C','A','G','A','C','T','G','C','G','C','G','C','A','A','G','T','C','A','A','C','G','G','G','T','G','G','C','A','A','G','G','C','G','C','C','G','C','G','C','A','A','G','C','A','G','C','T','G','G','C','C','A','C','C','A','A','G','G','T','G','G','C','T','C','G','C','A','A','G','A','G','C','G','C','A','C','C','T','G','C','C','A','C','T','G','G','C','G','G','C','G','T','G','A','A','G','A','A','G','C','C','G','C','A','C','C','G','C','T','A','C','C','G','G','C','C','C','G','G','C','A','C','G','G','T','G','G','C','G','C','T','T','C','G','C','G','A','G','A','T','C','C','G','C','C','G','C','T','A','C','C','A','G','A','A','G','T','C','C','A','C','T','G','A','G','C','T','G','C','T','A','A','T','C','C','G','C','A','A','G','T','T','G','C','C','C','T','T','C','C','A','G','C','G','G','C','T','G','A','T','G','C','G','C','G','A','G','A','T','C','G','C','T','C','A','G','G','A','C','T','T','T','A','A','G','A','C','C','G','A','C','C','T','G','C','G','C','T','T','C','C','A','G','A','G','C','T','C','G','G','C','C','G','T','G','A','T','G','G','C','G','C','T','G','C','A','G','G','A','G','G','C','G','T','G','C','G','A','G','T','C','T','T','A','C','C','T','G','G','T','G','G','G','G','C','T','G','T','T','T','G','A','G','G','A','C','A','C','C','A','A','C','C','T','G','T','G','T','G','T','C','A','T','C','C','A','T','G','C','C','A','A','A','C','G','G','G','T','C','A','C','C','A','T','C','A','T','G','C','C','T','A','A','G','G','A','C','A','T','C','C','A','G','C','T','G','G','C','A','C','G','C','C','G','T','A','T','C','C','G','C','G','G','G','G','A','G','C','G','G','G','C','C','T','A','G','G','A','G','G','G','C','T','A','T','C','T','C','G','C','C','A','C','C','T','G','A','G','A','G','G','T','T','G','C','G','C','A','A','C','G','T','T','C','A','C','C','C','C','A','A','A','G','G','C','T','C','T','T','T','T','A','A','G','A','G','C','C','A','C','C','C','A','C','C','T')

start_probability = {'exon': 0.50, 'donor1': 0.00, 'donor2': 0.00, 'intron': 0.50, 'acceptor1': 0.00, 'acceptor2': 0.00}
 
transition_probability = {
   'exon'      : {'exon': 0.90, 'donor1': 0.10, 'donor2': 0.00, 'intron': 0.00, 'acceptor1': 0.00, 'acceptor2': 0.00},
   'donor1'    : {'exon': 0.00, 'donor1': 0.00, 'donor2': 1.00, 'intron': 0.00, 'acceptor1': 0.00, 'acceptor2': 0.00},
   'donor2'    : {'exon': 0.00, 'donor1': 0.00, 'donor2': 0.00, 'intron': 1.00, 'acceptor1': 0.00, 'acceptor2': 0.00},
   'intron'    : {'exon': 0.00, 'donor1': 0.00, 'donor2': 0.00, 'intron': 0.90, 'acceptor1': 0.10, 'acceptor2': 0.00},
   'acceptor1' : {'exon': 0.00, 'donor1': 0.00, 'donor2': 0.00, 'intron': 0.00, 'acceptor1': 0.00, 'acceptor2': 1.00},
   'acceptor2' : {'exon': 1.00, 'donor1': 0.00, 'donor2': 0.00, 'intron': 0.00, 'acceptor1': 0.00, 'acceptor2': 0.00}
   }
 
emission_probability = {
   'exon'      : {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
   'donor1'    : {'A': 0.05, 'C': 0.00, 'G': 0.95, 'T': 0.00},
   'donor2'    : {'A': 0.00, 'C': 0.05, 'G': 0.00, 'T': 0.95},
   'intron'    : {'A': 0.40, 'C': 0.10, 'G': 0.10, 'T': 0.40},
   'acceptor1' : {'A': 0.95, 'C': 0.00, 'G': 0.05, 'T': 0.00},
   'acceptor2' : {'A': 0.05, 'C': 0.00, 'G': 0.95, 'T': 0.00}
   }
'''

states = ('A+', 'T+', 'C+', 'G+', 'A-', 'T-', 'G-', 'C-')
 
observations = []
with open('chr_22.txt') as f:
  while True:
    c = f.read(1)
    if not c:
      break 
    observations.append(c)
 
start_probability = {'A+': 0.125, 'T+': 0.125, 'G+': 0.125, 'C+': 0.125, 'A-': 0.125, 'T-': 0.125, 'G-': 0.125, 'C-': 0.125}
 
###################### VALUE TO CHANGE ############################    
##################### TRANSITION WEIGHT ########################### 
p = 0.1


transition_probability = {
   'A+' : {'A+': 0.180 *  (1-p), 'C+': 0.268 *  (1-p), 'G+': 0.430 *  (1-p), 'T+': 0.122 *  (1-p), 'A-': 0.25 * p, 'C-': 0.25 * p, 'G-': 0.25 * p, 'T-': 0.25 * p},
   'C+' : {'A+': 0.191 *  (1-p), 'C+': 0.299 *  (1-p), 'G+': 0.299 *  (1-p), 'T+': 0.211 *  (1-p), 'A-': 0.25 * p, 'C-': 0.25 * p, 'G-': 0.25 * p, 'T-': 0.25 * p},
   'G+' : {'A+': 0.161 *  (1-p), 'C+': 0.346 *  (1-p), 'G+': 0.373 *  (1-p), 'T+': 0.120 *  (1-p), 'A-': 0.25 * p, 'C-': 0.25 * p, 'G-': 0.25 * p, 'T-': 0.25 * p},
   'T+' : {'A+': 0.082 *  (1-p), 'C+': 0.357 *  (1-p), 'G+': 0.391 *  (1-p), 'T+': 0.170 *  (1-p), 'A-': 0.25 * p, 'C-': 0.25 * p, 'G-': 0.25 * p, 'T-': 0.25 * p},
   

   'A-' : {'A+': 0.25 * p, 'C+': 0.25 * p, 'G+': 0.25 * p, 'T+': 0.25 * p, 'A-': 0.300 *  (1-p), 'C-': 0.200 *  (1-p), 'G-': 0.290 *  (1-p), 'T-': 0.210 *  (1-p)},
   'C-' : {'A+': 0.25 * p, 'C+': 0.25 * p, 'G+': 0.25 * p, 'T+': 0.25 * p, 'A-': 0.319 *  (1-p), 'C-': 0.302 *  (1-p), 'G-': 0.081 *  (1-p), 'T-': 0.291 *  (1-p)},   
   'G-' : {'A+': 0.25 * p, 'C+': 0.25 * p, 'G+': 0.25 * p, 'T+': 0.25 * p, 'A-': 0.251 *  (1-p), 'C-': 0.251 *  (1-p), 'G-': 0.299 *  (1-p), 'T-': 0.199 *  (1-p)},
   'T-' : {'A+': 0.25 * p, 'C+': 0.25 * p, 'G+': 0.25 * p, 'T+': 0.25 * p, 'A-': 0.176 *  (1-p), 'C-': 0.242 *  (1-p), 'G-': 0.291 *  (1-p), 'T-': 0.291 *  (1-p)}

   }
 
emission_probability = {
   'A+' : {'A': 1.0, 'T': 0.0, 'G': 0.0, 'C': 0.0},
   'T+' : {'A': 0.0, 'T': 1.0, 'G': 0.0, 'C': 0.0},
   'G+' : {'A': 0.0, 'T': 0.0, 'G': 1.0, 'C': 0.0},
   'C+' : {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 1.0},
   'A-' : {'A': 1.0, 'T': 0.0, 'G': 0.0, 'C': 0.0},
   'T-' : {'A': 0.0, 'T': 1.0, 'G': 0.0, 'C': 0.0},
   'G-' : {'A': 0.0, 'T': 0.0, 'G': 1.0, 'C': 0.0},
   'C-' : {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 1.0}
   }
   
###################################################
############  END MODIFICATIONS HERE  #############
###################################################  

def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    path = {}
 
    # Initialize base cases (t == 0)
    for y in states:
        V[0][y] = start_p[y] * emit_p[y][obs[0]]
        path[y] = [y]
 
    # Run Viterbi for t > 0
    for t in range(1, len(obs)):
        V.append({})
        newpath = {}
 
        for y in states:
            (prob, state) = max((V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states)
            V[t][y] = prob
            newpath[y] = path[state] + [y]
 
        # Don't need to remember the old paths
        path = newpath
    n = 0           # if only one element is observed max is sought in the initialization values
    if len(obs) != 1:
        n = t
    print_dptable(V)
    (prob, state) = max((V[n][y], y) for y in states)
    return (prob, path[state])
 
# Don't study this, it just prints a table of the steps.
def print_dptable(V):
    s = "    " + " ".join(("%7d" % i) for i in range(len(V))) + "\n"
    for y in V[0]:
        s += "%.5s: " % y
        s += " ".join("%.7s" % ("%f" % v[y]) for v in V)
        s += "\n"
    print(s)

def example():
    return viterbi(observations,
                   states,
                   start_probability,
                   transition_probability,
                   emission_probability)
print(example())