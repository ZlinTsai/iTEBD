# Infinite Time Evolving Block Decimation (iTEBD)

## Name: Zheng-Lin Tsai             &ensp; Student ID: s106022501

### Method

In quantum mechanics, the time-dependent wave function is <img src="https://latex.codecogs.com/png.latex?\left&space;|&space;\Psi(t)&space;\right&space;\rangle&space;=&space;U(t)\left&space;|&space;\Psi(0)&space;\right&space;\rangle" title="\left | \Psi(t) \right \rangle = U(t)\left | \Psi(0) \right \rangle" /> , U(t) is time evolution operator. <br />

By Schrödinger equation, we can get  <img src="https://latex.codecogs.com/png.latex?\large&space;U(t)&space;=&space;e^{-i\frac{\hat{H}}{\hbar}t}" title="\large U(t) = e^{-i\frac{\hat{H}}{\hbar}t}" />, then set <img src="https://latex.codecogs.com/png.latex?\hbar" title="\hbar" />=1 and use imaginary time <img src="https://latex.codecogs.com/png.latex?\tau" title="\tau" />, so rewrite imaginary time operator <img src="https://latex.codecogs.com/gif.latex?\large&space;U(\tau)&space;=&space;e^{-\hat{H}\tau}" title="\large U(\tau) = e^{-\hat{H}\tau}" /> . <br />

Next, <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}&space;=&space;h^{[1,2]}&space;&plus;&space;h^{[2,3]}&space;&plus;&space;h^{[3,4]}&space;&plus;&space;..." title="\large \hat{H} = h^{[1,2]} + h^{[2,3]} + h^{[3,4]} + ..." /> h is local two site hamiltonian and [1,2] is index of site , collect [1,2], [3,4], [5,6]... into odd hamiltonian <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}_{odd}" title="\large \hat{H}_{odd}" /> , also <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}_{even}" title="\large \hat{H}_{even}" />. <br />
So <img src="https://latex.codecogs.com/gif.latex?\large&space;\hat{H}&space;=&space;\hat{H}_{odd}&space;&plus;&space;\hat{H}_{even}" title="\large \hat{H} = \hat{H}_{odd} + \hat{H}_{even}" />.  <br />

And use Trotter-Suzuki decomposition <img src="https://latex.codecogs.com/gif.latex?\large&space;e^{(A&plus;B)\delta}&space;\approx&space;e^{A\delta}e^{B\delta}&space;&plus;&space;O(\delta^2)" title="\large e^{(A+B)\delta} \approx e^{A\delta}e^{B\delta} + O(\delta^2)" /> , <img src="https://latex.codecogs.com/gif.latex?\large&space;\delta" title="\large \delta" /> is a small parameter. <br />

In order to use Trotter-Suzuki decomposition, we need slice imaginary time <img src="https://latex.codecogs.com/gif.latex?\large&space;\tau&space;=&space;\frac{\tau}{N}N&space;\equiv&space;\delta&space;N" title="\large \tau = \frac{\tau}{N}N \equiv \delta N" /> , N is number of slices. <br />

So we get <img src="https://latex.codecogs.com/gif.latex?\large&space;e^{-\hat{H}\tau}&space;=&space;\prod^{N}&space;e^{-\hat{H}\delta}" title="\large e^{-\hat{H}\tau} = \prod^{N} e^{-\hat{H}\delta}" /> , and <img src="https://latex.codecogs.com/gif.latex?\large&space;e^{-\hat{H}\delta}&space;\approx&space;e^{-\hat{H}_{odd}\delta}e^{-\hat{H}_{even}\delta}" title="\large e^{-\hat{H}\delta} \approx e^{-\hat{H}_{odd}\delta}e^{-\hat{H}_{even}\delta}" /> . In ![FIG.1](database/fig_1.png)
