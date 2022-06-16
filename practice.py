import numpy as np
class Solution:
    def CRR_method(self, K,T,S0,r,N,sigma,opttype='C'):
    #precomute constants
        dt = T/N
        u = np.exp(sigma*np.sqrt(dt))
        d = 1/u
        q = (np.exp(r*dt) - d) / (u-d)
        # print(q)
        disc = np.exp(-r*dt)
        
        # initialise asset prices at maturity - Time step N
        S = np.zeros(N+1)
        S[0] = S0*d**N
        # print(S[0])
        for j in range(1,N+1):
            S[j] = S[j-1]*u/d
            
        # print(S)
        
        # initialise option values at maturity
        C = np.zeros(N+1)
        for j in range(0,N+1):
            if opttype == 'C':
                C[j] = max(0, S[j]-K)
            else:
                C[j] = max(0, K - S[j])
                
        # print(C)
            
        # step backwards through tree
        for i in np.arange(N,0,-1):
            for j in range(0,i):
                # print(C[j])
                C[j] = disc * ( q*C[j+1] + (1-q)*C[j] )
                # print(C[j])
                # print('')
            # print(C)
        return C[0]

    def american_slow_tree(self, K,T,S0,r,N,sigma,opttype='P'):
        #precompute values
        dt = T/N
        u = np.exp(sigma*np.sqrt(dt))
        d = 1/u
        q = (np.exp(r*dt) - d)/(u-d)
        disc = np.exp(-r*dt)
        
        # initialise stock prices at maturity
        S = np.zeros(N+1)
        for j in range(0, N+1):
            S[j] = S0 * u**j * d**(N-j)
            
        # option payoff 
        C = np.zeros(N+1)
        for j in range(0, N+1):
            if opttype == 'P':
                C[j] = max(0, K - S[j])
            else:
                C[j] = max(0, S[j] - K)
        
        # backward recursion through the tree
        for i in np.arange(N-1,-1,-1):
            for j in range(0,i+1):
                S = S0 * u**j * d**(i-j)
                C[j] = disc * ( q*C[j+1] + (1-q)*C[j] )
                if opttype == 'P':
                    C[j] = max(C[j], K - S)
                else:
                    C[j] = max(C[j], S - K)
                    
        return C[0]

s = Solution();
K = 100;
T = 1;
S0 = 100;
r = 0.05;
N = 10;
sigma = 0.2;

print(s.CRR_method(K,T,S0,r,N,sigma,opttype='C'))
print()
print(s.american_slow_tree(K,T,S0,r,N,sigma,opttype='C'))