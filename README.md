# TruenorthSNN-Simulator
This is a simulator for IBM Truenorth SNN(Spike Neuron Network) architecture implemented in MATLAB code.The neuron specification equations in this code are based on an IBM article *"Cognitive Computing Building Block: A Versatile and Efficient Digital Neuron Model for Neurosynaptic Cores"*.   
## Manual
1. Download this repository as a folder and open it in MATLAB or Octave environment.  
  
2. Type `top` in the command line. It should show this:  
  ```
  Truenorth SNN Architecture Simulator V0.02 
  [1] Configure a new network
  [2] Use an old configuration file(.mat)
  ```
  Input 1 or 2 to select a mode, and then follow the instructions printed in the command line.  
  
  **Notice 1:** `Time interval` can be set to from 0 to 1000. The bigger you make this variable, the more accurate your result will be, the more time the calculation will spend.  
  
  **Notice 2:** You should input data in MATLAB style. For example, if you want to input a vector, you should surround the numbers with `[]` and seperate them with a space or a comma, like `[1 2]`.  
  
  **Notice 3:** The following prompt prompts you to input the spike rate expressions inputed by axons. The `x` in the expression will be a `1:max_rate` vector in the code. As `max_rate` is set to 1000 by default, it is recommanded to keep the expression's max value below 1000. Spike rate above 1000Hz will have the same effect as 1000Hz.   
  ```
  Input the coefficient vector of spike rate expression for every axon.
  Linear expressions only. For example, if the expression is x+2, input [1 2].
  Spike rate expression coefficient vector of axon 1:
  ```
  You'd better input something like `[0.7 5]`.  
  
3. When the calculation is finished, you can save the current variable table as a `.mat` file when the following prompt appears, in order to save your neuron network configuration.  
  `Save current configuration file?[y/n]`  
  And next time when you start this simulator, you can reuse this configuration file.  
  
4. There are a few configuration files in the `./omat` and `./mmat` folder which can perform arithmetic functions. You can test them by choosing the second mode and then type `./omat/(mat file name)` or `./mmat/(mat file name)`.  

5. The simulator will plot a figure indicating the input spike rate and output spike rate. The input curve is in green and the output curve is in red by default. The `plot` function receives 3 vectors(2 inputs and 1 output) in this code.  
  **Notice 4:** You may have to manually modify the `plot` function parameters to fit your own input and output.  
  
## Configuration File Format
The configuration file is simply a `.mat` file, generated by `save` command, but there is some difference between MATLAB's mat file and Octave's mat file. Both file saves current workspace's variable table. However, MATAB's mat file is an encrypted binary file which can be read by both MATLAB and Octave, and Octave's mat file is a plain text file but can only be read by Octave. For the convenience of mat file modification, it is recommanded to use Octave and its mat file.  

MATLAB's mat files are saved in `./mmat`, and Octave's mat files are saved in `./omat`.  

## Example
This example will generate a SNN that can preform arithmetic addition.
```
Truenorth SNN Architecture Simulator V0.02
[1] Configure a new network
[2] Use an old configuration file(.mat)
1
Axon number :
2
Neuron number :
1
Time interval :
500
Synapse connection vector of neuron 1 with 2 axons:
[1 1]
Synapse mode vector of neuron 1's 2 synapses:
[0 0]
Synapse weight vector of neuron 1's 2 synapses:
[1 1]
Leak mode vector of 1 neurons:
0
Leak weight vector of 1 neurons:
0
Leak-reversal flag vector of 1 neurons:
0
Reset mode vector of 1 neurons:
1
Reset voltage vector of 1 neurons:
0
Positive threshold vector of 1 neurons:
1
Negative threshold vector of 1 neurons:
1
Threshold PRN mask vector of 1 neurons(0 to 255):
0
Negative thresh mode vector of 1 neurons:
0
Input the coefficient vector of spike rate expression for every axon.
Linear expressions only. For example, if the expression is x+2, input [1 2]. 
Spike rate expression coefficient vector of axon 1:
[0.2 1]
Spike rate expression coefficient vector of axon 2:
[0.5 3]
```

Or use an old configuration file.  
```
Truenorth SNN Architecture Simulator V0.02
[1] Configure a new network
[2] Use an old configuration file(.mat)
2
Input the configuration file name:
./mmat/add.mat
Input the coefficient vector of spike rate expression for every axon.
Linear expressions only. For example, if the expression is x+2, input [1 2]. 
Spike rate expression coefficient vector of axon 1:
[0.2 1]
Spike rate expression coefficient vector of axon 2:
[0.5 3]

```
Output figure of this example:  
![Result](http://p1.bpimg.com/1949/09d7e023df88d0ce.png)
