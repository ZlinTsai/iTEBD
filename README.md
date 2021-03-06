# Infinite Time Evolving Block Decimation (iTEBD)

## Name: Zheng-Lin Tsai             &ensp; Student ID: s106022501

### 1.Method

In quantum mechanics, the time-dependent wave function is <img src="https://latex.codecogs.com/png.latex?\left&space;|&space;\Psi(t)&space;\right&space;\rangle&space;=&space;U(t)\left&space;|&space;\Psi(0)&space;\right&space;\rangle" title="\left | \Psi(t) \right \rangle = U(t)\left | \Psi(0) \right \rangle" /> , U(t) is time evolution operator. <br />

If time is very long, we expect state is most stable state. So we can find ground state <img src="https://latex.codecogs.com/gif.latex?\large&space;\left&space;|&space;\Psi_{GS}&space;\right&space;\rangle&space;=&space;\lim_{t&space;\rightarrow&space;\infty&space;}&space;e^{-i\frac{\hat{H}}{\hbar}t}&space;\left&space;|&space;\Psi(0)&space;\right&space;\rangle" title="\large \left | \Psi_{GS} \right \rangle = \lim_{t \rightarrow \infty } e^{-i\frac{\hat{H}}{\hbar}t} \left | \Psi(0) \right \rangle" /> . <br />

By Schrödinger equation, we can get  <img src="https://latex.codecogs.com/png.latex?\large&space;U(t)&space;=&space;e^{-i\frac{\hat{H}}{\hbar}t}" title="\large U(t) = e^{-i\frac{\hat{H}}{\hbar}t}" />, then set <img src="https://latex.codecogs.com/png.latex?\hbar" title="\hbar" />=1 and use imaginary time <img src="https://latex.codecogs.com/png.latex?\tau" title="\tau" />, so rewrite imaginary time operator <img src="https://latex.codecogs.com/gif.latex?\large&space;U(\tau)&space;=&space;e^{-\hat{H}\tau}" title="\large U(\tau) = e^{-\hat{H}\tau}" /> . <br />

Next, <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}&space;=&space;h^{(1,2)}&space;&plus;&space;h^{(2,3)}&space;&plus;&space;h^{(3,4)}&space;&plus;&space;..." title="\large \hat{H} = h^{[1,2]} + h^{[2,3]} + h^{[3,4]} + ..." /> h is local two site hamiltonian and (1,2) is index of site , collect (1,2), (3,4), (5,6)... into odd hamiltonian <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}_{odd}" title="\large \hat{H}_{odd}" /> , also <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}_{even}" title="\large \hat{H}_{even}" />. <br />
So <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}&space;=&space;\hat{H}_{odd}&space;&plus;&space;\hat{H}_{even}" title="\large \hat{H} = \hat{H}_{odd} + \hat{H}_{even}" /> [1, 2].  <br />

And use Trotter-Suzuki decomposition <img src="https://latex.codecogs.com/gif.latex?\large&space;e^{(A&plus;B)\delta}&space;\approx&space;e^{A\delta}e^{B\delta}&space;&plus;&space;O(\delta^2)" title="\large e^{(A+B)\delta} \approx e^{A\delta}e^{B\delta} + O(\delta^2)" /> , <img src="https://latex.codecogs.com/gif.latex?\large&space;\delta" title="\large \delta" /> is a small parameter. <br />

In order to use Trotter-Suzuki decomposition, we need slice imaginary time <img src="https://latex.codecogs.com/gif.latex?\large&space;\tau&space;=&space;\frac{\tau}{N}N&space;\equiv&space;\delta&space;N" title="\large \tau = \frac{\tau}{N}N \equiv \delta N" /> , N is number of slices. <br />

So we get <img src="https://latex.codecogs.com/gif.latex?\large&space;e^{-\hat{H}\tau}&space;=&space;\prod^{N}&space;e^{-\hat{H}\delta}" title="\large e^{-\hat{H}\tau} = \prod^{N} e^{-\hat{H}\delta}" /> , and <img src="https://latex.codecogs.com/gif.latex?\large&space;e^{-\hat{H}\delta}&space;\approx&space;e^{-\hat{H}_{odd}\delta}e^{-\hat{H}_{even}\delta}" title="\large e^{-\hat{H}\delta} \approx e^{-\hat{H}_{odd}\delta}e^{-\hat{H}_{even}\delta}" /> , show graphically ![FIG.1](database/fig_1.PNG)  <br />

Above picture express by Matrix Product States (MPS)[1], and we consider infinite chain that mean each site must the same (translation invariant), we only updata two site and then we can get ground state. <br />

### 2.Encoding
The code [iTEBD_Ising.cpp](iTEBD_Ising.cpp) use Uni10 library[3]. <br />
![FIG.2](database/fig_2.PNG) <br />
The iTEBD update scheme, <br />
(i) contract initial MPS into <img src="https://latex.codecogs.com/gif.latex?\large&space;\Theta" title="\large \Theta" />, <br />
(ii) applying U = <img src="https://latex.codecogs.com/gif.latex?\large&space;e^{-\hat{H}\delta}" title="\large e^{-\hat{H}\delta}" />, <br />
(iii) SVD <img src="https://latex.codecogs.com/gif.latex?\large&space;\widetilde{\Theta}" title="\large \widetilde{\Theta}" /> in order to get the same form with (i), <br />
(iv) apply origin <img src="https://latex.codecogs.com/gif.latex?\large&space;\Lambda^B" title="\large \Lambda^B" /> and inverse <img src="https://latex.codecogs.com/gif.latex?\large&space;(\Lambda^B)^{-1}" title="\large (\Lambda^B)^{-1}" /> , then contract <img src="https://latex.codecogs.com/gif.latex?\large&space;(\Lambda^B)^{-1}" title="\large (\Lambda^B)^{-1}" /> with X and Y, so get origin form after updating. <br />

And also even term (exchange A and B site), repeat (i)~(iv) until energy convergence. <br />

### 3.Application
Consider one-Dimension Ising model in tranverse magnetic field  <br />

<img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}&space;=&space;J\sum_i&space;\sigma^z_i&space;\sigma^z_{i&plus;1}&space;&plus;&space;h\sum_i&space;\sigma^x_i" title="\large \hat{H} = J\sum_i \sigma^z_i \sigma^z_{i+1} + h\sum_i \sigma^x_i" />  <br /> 


 <img src="https://latex.codecogs.com/gif.latex?\large&space;\sigma" title="\large \sigma" />  is Pauli matrix. <br /> 
 
![FIG.2](database/fig_3.PNG) <br />
 
When h field is small system is antiferromagnetic in z-axis by J > 0, if J < 0 system will be ferromagnetic. When h field is large, system is paramagnetic. And h = 1 is critical point.
 
![FIG.4](database/fig_4.PNG) <br />
 
![FIG.5](database/fig_5.PNG) <br />
 
![FIG.6](database/fig_6.PNG) <br />

### 4.References
[1] Frank Pollmann Efficient Numerical Simulations Using Matrix-Product States  <br />
[2] G. Vidal. Classical simulation of infinite-size quantum lattice systems in one spatial dimension. arXiv:cond-mat/0605597   <br />
[3] Uni10: an open-source library for tensor network algorithms. https://gitlab.com/uni10/uni10  <br />
