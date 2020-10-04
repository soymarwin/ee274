%% Marwin B. Alejo   2020-20221   EE274_ProgEx02
% 

%% System 1 : $y[n] = 0.5x[n] + 0.5x[n-1]$
% Using recursive method.

[x1,fs_x1] = audioread('inputs/x1.wav'); 
x1_1 = x1(:,1); 
x1_2 = x1(:,2);
[x2,fs_x2] = audioread('inputs/x2.wav');
[x3,fs_x3] = audioread('inputs/x3.wav');
[x4,fs_x4] = audioread('inputs/x4.wav');
[x5,fs_x5] = audioread('inputs/x5.wav');

% Input 1 channel 1 to System 1
    y1_1 = dt_1(x1_1);
    figure();
    subplot 211
    stem(1:length(x1_1),x1_1); title('input signal');
    subplot 212
    stem(1:length(y1_1),y1_1); title('output signal');

% Input 1 channel 2 to System 1
    y1_2 = dt_1(x1_2);
    figure();
    subplot 211
    stem(1:length(x1_2),x1_2); title('input signal');
    subplot 212
    stem(1:length(y1_2),y1_2); title('output signal');
y1_1_both = [y1_1(:),y1_2(:)]; %combine two channels.
% soundsc(x1,fs_x1); %x1 original audio
% soundsc(y1_1_both,fs_x1); %x1 output audio through system 1

%%
% _OBSERVATION: Although there is a clear small difference between the input
% and output of audio x1 as per the system 1 generated graphs, it cannot be
% sense through bare hearing._
    
% Input 2 to System 1
y1_3 = dt_1(x2);
figure();
subplot 211
stem(1:length(x2),x2); title('input signal');
subplot 212
stem(1:length(y1_3),y1_3); title('output signal');
% soundsc(x2,fs_x2); %x2 original audio
% soundsc(y1_3,fs_x2); %x2 output audio through system 1

%%
% _OBSERVATION: Although there is a clear small difference between the input
% and output of audio x2 as per the system 1 generated graphs, it cannot be
% sense through bare hearing._

% Input 3 to System 1
y1_4 = dt_1(x3);
figure();
subplot 211
stem(1:length(x3),x3); title('input signal');
subplot 212
stem(1:length(y1_4),y1_4); title('output signal');
% soundsc(x3,fs_x3); %x3 original audio
% soundsc(y1_4,fs_x3); %x3 output audio through system 1

%%
% _OBSERVATION: Although there is a clear small difference between the input
% and output of audio x3 as per the system 1 generated graphs, it cannot be
% sense through bare hearing._

% Input 4 to System 1
y1_5 = dt_1(x4);
figure();
subplot 211
stem(1:length(x4),x4); title('input signal');
subplot 212
stem(1:length(y1_5),y1_5); title('output signal');
% soundsc(x4,fs_x4); %x4 original audio
% soundsc(y1_5,fs_x4); %x4 output audio through system 1

%%
% _OBSERVATION: Although there is a clear small difference between the input
% and output of audio x4 as per the system 1 generated graphs, it cannot be
% sense through bare hearing._

% Input 5 to System 1
y1_6 = dt_1(x5);
figure();
subplot 211
stem(1:length(x5),x5); title('input signal');
subplot 212
stem(1:length(y1_6),y1_6); title('output signal');
% soundsc(x5,fs_x5); %x5 original audio
% soundsc(y1_6,fs_x5); %x5 output audio through system 1

%%
% _OBSERVATION: Although there is a clear small difference between the input
% and output of audio x5 as per the system 1 generated graphs, it cannot be
% sense through bare hearing._

%%
% _OVERALL, the input and output signals through system 1 greatly oes not
% differ from each other. The signals are both audible in both ends._

%% System 2 : $y[n] = x[n] - 2y[n-1] - 2y[n-2]$
% Using recursive method.

% Input 1 channel 1 to System 2
    y2_1 = dt_2(x1_1);
    figure();
    subplot 211
    stem(1:length(x1_1),x1_1); title('input signal');
    subplot 212
    stem(1:length(y2_1),y2_1); title('output signal');
% Input 1 channel 2 to System 2
    y2_2 = dt_2(x1_2);
    figure();
    subplot 211
    stem(1:length(x1_2),x1_2); title('input signal');
    subplot 212
    stem(1:length(y2_2),y2_2); title('output signal');
y2_1_both = [y2_1(:),y2_2(:)]; %combine two channels.
% soundsc(x1,fs_x1); %x1 original audio
% soundsc(y2_1_both,fs_x1); %x1 output audio through system 2

%%
% _OBSERVATION: As per the system 2 generated graph and auditory comparison
% of the input and output signals of x1, there exist a great difference
% between the two as such the amplitude of the output signal was compressed
% from 0 to near end and was compressed in time-domain at the end._

% Input 2 to System 2
y2_3 = dt_2(x2);
figure();
subplot 211
stem(1:length(x2),x2); title('input signal');
subplot 212
stem(1:length(y2_3),y2_3); title('output signal');
% soundsc(x2,fs_x2); %x2 original audio
% soundsc(y2_3,fs_x2); %x2 output audio through system 2

%%
% _OBSERVATION: As per the system 2 generated graph and auditory comparison
% of the input and output signals of x2, there exist a great difference
% between the two as such the amplitude of the output signal was compressed
% from 0 to near end and was compressed in time-domain at the end._

% Input 3 to System 2
y2_4 = dt_2(x3);
figure();
subplot 211
stem(1:length(x3),x3); title('input signal');
subplot 212
stem(1:length(y2_4),y2_4); title('output signal');
% soundsc(x3,fs_x3); %x3 original audio
% soundsc(y2_4,fs_x3); %x3 output audio through system 2

%%
% _OBSERVATION: As per the system 2 generated graph and auditory comparison
% of the input and output signals of x3, there exist a great difference
% between the two as such the amplitude of the output signal was compressed
% from 0 to near end and was compressed in time-domain at the end._

% Input 4 to System 2
y2_5 = dt_2(x4);
figure();
subplot 211
stem(1:length(x4),x4); title('input signal');
subplot 212
stem(1:length(y2_5),y2_5); title('output signal');
% soundsc(x4,fs_x4); %x4 original audio
% soundsc(y2_5,fs_x4); %x4 output audio through system 2

%%
% _OBSERVATION: As per the system 2 generated graph and auditory comparison
% of the input and output signals of x4, there exist a great difference
% between the two as such the amplitude of the output signal was compressed
% from 0 to near end and was compressed in time-domain at the end._

% Input 5 to System 2
y2_6 = dt_2(x5);
figure();
subplot 211
stem(1:length(x5),x5); title('input signal');
subplot 212
stem(1:length(y2_6),y2_6); title('output signal');
% soundsc(x5,fs_x5); %x5 original audio
% soundsc(y2_6,fs_x5); %x5 output audio through system 2

%%
% _OBSERVATION: As per the system 2 generated graph and auditory comparison
% of the input and output signals of x5, there exist a great difference
% between the two as such the amplitude of the output signal was compressed
% from 0 to near end and was compressed in time-domain at the end._

%%
% _OVERALL, the input signals through system 2 were compressed in terms of
% space-domain resulting for output signals to expand in time-domain. The
% signals are no longer audible on human-cognitive level._

%% System 3 : $y[n] = 1.5x[n] + 0.5x[n-1] -2y[n-1] - 2y[n-2]$
% Using impulse response and convolution.

% Input 1 channel 1 to System 3
    y3_1 = dt_3(x1_1);
    figure();
    subplot 211
    stem(1:length(x1_1),x1_1); title('input signal');
    subplot 212
    stem(1:length(y3_1),y3_1); title('output signal');
% Input 1 channel 2 to System 3
    y3_2 = dt_3(x1_2);
    figure();
    subplot 211
    stem(1:length(x1_2),x1_2); title('input signal');
    subplot 212
    stem(1:length(y3_2),y3_2); title('output signal');
y3_1_both = [y3_1(:),y3_2(:)]; %combine two channels.
% soundsc(x1,fs_x1); %x1 original audio
% soundsc(y3_1_both,fs_x1); %x1 output audio through system 3

%%
% _OBSERVATION: The output signal of x1 is smoother than its input after
% passing system 3._

% Input 2 to System 3
y3_3 = dt_3(x2);
figure();
subplot 211
stem(1:length(x2),x2); title('input signal');
subplot 212
stem(1:length(y3_3),y3_3); title('output signal');
% soundsc(x2,fs_x2); %x2 original audio
% soundsc(y3_3,fs_x2); %x2 output audio through system 3

%%
% _OBSERVATION: The output signal of x2 is smoother than its input after
% passing system 3._

% Input 3 to System 3
y3_4 = dt_3(x3);
figure();
subplot 211
stem(1:length(x3),x3); title('input signal');
subplot 212
stem(1:length(y3_4),y3_4); title('output signal');
% soundsc(x3,fs_x3); %x3 original audio
% soundsc(y3_4,fs_x3); %x3 output audio through system 3

%%
% _OBSERVATION: The output signal of x3 is smoother than its input after
% passing system 3._

% Input 4 to System 3
y3_5 = dt_3(x4);
figure();
subplot 211
stem(1:length(x4),x4); title('input signal');
subplot 212
stem(1:length(y3_5),y3_5); title('output signal');
% soundsc(x4,fs_x4); %x4 original audio
% soundsc(y3_5,fs_x4); %x4 output audio through system 3

%%
% _OBSERVATION: The output signal of x4 is smoother than its input after
% passing system 3._

% Input 5 to System 3
y3_6 = dt_3(x5);
figure();
subplot 211
stem(1:length(x5),x5); title('input signal');
subplot 212
stem(1:length(y3_6),y3_6); title('output signal');
% soundsc(x5,fs_x5); %x5 original audio
% soundsc(y3_6,fs_x5); %x5 output audio through system 3

%%
% _OBSERVATION: The output signal of x5 is smoother than its input after
% passing system 3._

%%
% _OVERALL, output signals after passing system 3 become smoother in audio
% as compared with its input counterpart. The output signals amplitude, as
% observed through charts are smaller by a small percentage.__

%% System 4 : $y[n] = x[n] + 0.5y[n-L] + 0.5y[n-L-1]$
% Using impulse response and convolution.

% Input 1 channel 1 to System 4
    y4_1_50 = dt_4(x1_1,50); % x1 L=50 ch1 system4
    y4_1_100 = dt_4(x1_1,100); % x1 L=100 ch1 system 4
    y4_1_400 = dt_4(x1_1,400); % x1 L=400 ch1 system 4
    figure();
    subplot 411
    stem(1:length(x1_1),x1_1); title('input signal');
    subplot 412
    stem(1:length(y4_1_50),y4_1_50); title('x1 ch1 L=50 output signal');
    subplot 413
    stem(1:length(y4_1_100),y4_1_100); title('x1 ch1 L=100 output signal');
    subplot 414
    stem(1:length(y4_1_400),y4_1_400); title('x1 ch1 L=400 output signal');
% Input 1 channel 2 to System 4
    y4_2_50 = dt_4(x1_2,50); % x1 L=50 ch2 system4
    y4_2_100 = dt_4(x1_2,100); % x1 L=100 ch2 system4
    y4_2_400 = dt_4(x1_2,400); % x1 L=400 ch2 system4
    figure();
    subplot 411
    stem(1:length(x1_2),x1_2); title('input signal');
    subplot 412
    stem(1:length(y4_2_50),y4_2_50); title('x1 ch2 L=50 output signal');
    subplot 413
    stem(1:length(y4_2_100),y4_2_100); title('x1 ch2 L=100 output signal');
    subplot 414
    stem(1:length(y4_2_400),y4_2_400); title('x1 ch2 L=400 output signal');
y4_1_both_50 = [y4_1_50(:),y4_2_50(:)]; %combine two channels L=50
y4_1_both_100 = [y4_1_100(:),y4_2_100(:)]; %combine two channels L=100
y4_1_both_400 = [y4_1_400(:),y4_2_400(:)]; %combine two channels L=400
% soundsc(x1,fs_x1); %x1 original audio
% soundsc(y4_1_both_50,fs_x1); %x1 output audio through system 4 L=50
% soundsc(y4_1_both_100,fs_x1); %x1 output audio through system 4 L=100
% soundsc(y4_1_both_400,fs_x1); %x1 output audio through system 4 L=400

%%
% _OBSERVATION: Output signals after passing system 4 appear to be
% downsampled however, the output signal with L=4 has a clearer audio than
% those x1 with L=100 and L=50. Visually, input and output signals are hard
% to differ due to small percentage of difference._

% Input 2 to System 4
y4_3_50 = dt_4(x2,50); % x2 L=50 system4
y4_3_100 = dt_4(x2,100); % x2 L=100 system4
y4_3_400 = dt_4(x2,400); % x2 L=400 system4
figure();
subplot 411
stem(1:length(x2),x2); title('input signal');
subplot 412
stem(1:length(y4_3_50),y4_3_50); title('L=50 output signal');
subplot 413
stem(1:length(y4_3_100),y4_3_100); title('L=100 output signal');
subplot 414
stem(1:length(y4_3_400),y4_3_400); title('L=400 output signal');
% soundsc(x2,fs_x2); %x2 original audio
% soundsc(y4_3_50,fs_x2); %x2 output audio through system 4 L=50
% soundsc(y4_3_100,fs_x2); %x2 output audio through system 4 L=100
% soundsc(y4_3_400,fs_x2); %x2 output audio through system 4 L=400

%%
% _OBSERVATION: Output signals after passing system 4 appear to be
% downsampled however, the output signal with L=4 has a clearer audio than
% those x1 with L=100 and L=50. Visually, input and output signals are hard
% to differ due to small percentage of difference._

% Input 3 to System 4
y4_4_50 = dt_4(x3,50); % x3 L=50 system4
y4_4_100 = dt_4(x3,100); % x3 L=100 system4
y4_4_400 = dt_4(x3,400); % x3 L=400 system4
figure();
subplot 411
stem(1:length(x3),x3); title('input signal');
subplot 412
stem(1:length(y4_4_50),y4_4_50); title('L=50 output signal');
subplot 413
stem(1:length(y4_4_100),y4_4_100); title('L=100 output signal');
subplot 414
stem(1:length(y4_4_400),y4_4_400); title('L=400 output signal');
% soundsc(x3,fs_x3); %x3 original audio
% soundsc(y4_4_50,fs_x3); %x3 output audio through system 4 L=50
% soundsc(y4_4_100,fs_x3); %x3 output audio through system 4 L=100
% soundsc(y4_4_400,fs_x3); %x3 output audio through system 4 L=400

%%
% _OBSERVATION: Output signals after passing system 4 appear to be
% downsampled however, the output signal with L=4 has a clearer audio than
% those x1 with L=100 and L=50. Visually, input and output signals are hard
% to differ due to small percentage of difference._

% Input 4 to System 3
y4_5_50 = dt_4(x4,50); % x4 L=50 system4
y4_5_100 = dt_4(x4,100); % x4 L=100 system4
y4_5_400 = dt_4(x4,400); % x4 L=400 system4
figure();
subplot 411
stem(1:length(x4),x4); title('input signal');
subplot 412
stem(1:length(y4_5_50),y4_5_50); title('L=50 output signal');
subplot 413
stem(1:length(y4_5_100),y4_5_100); title('L=100 output signal');
subplot 414
stem(1:length(y4_5_400),y4_5_400); title('L=400 output signal');
% soundsc(x4,fs_x4); %x4 original audio
% soundsc(y4_5_50,fs_x4); %x4 output audio through system 4 L=50
% soundsc(y4_5_100,fs_x4); %x4 output audio through system 4 L=100
% soundsc(y4_5_400,fs_x4); %x4 output audio through system 4 L=400

%%
% _OBSERVATION: Output signals after passing system 4 appear to be
% downsampled however, the output signal with L=4 has a clearer audio than
% those x1 with L=100 and L=50. Visually, input and output signals are hard
% to differ due to small percentage of difference._

% Input 5 to System 4
y4_6 = dt_4(x5,100);
y4_6_50 = dt_4(x5,50); % x5 L=50 system4
y4_6_100 = dt_4(x5,100); % x5 L=100 system4
y4_6_400 = dt_4(x5,400); % x5 L=400 system4
figure();
subplot 411
stem(1:length(x5),x5); title('input signal');
subplot 412
stem(1:length(y4_6_50),y4_6_50); title('L=50 output signal');
subplot 413
stem(1:length(y4_6_100),y4_6_100); title('L=100 output signal');
subplot 414
stem(1:length(y4_6_400),y4_6_400); title('L=400 output signal');
% soundsc(x5,fs_x5); %x5 original audio
% soundsc(y4_6_50,fs_x5); %x5 output audio through system 4 L=50
% soundsc(y4_6_100,fs_x5); %x5 output audio through system 4 L=100
% soundsc(y4_6_400,fs_x5); %x5 output audio through system 4 L=400

%%
% _OBSERVATION: Output signals after passing system 4 appear to be
% downsampled however, the output signal with L=4 has a clearer audio than
% those x1 with L=100 and L=50. Visually, input and output signals are hard
% to differ due to small percentage of difference._

%%
% _OVERALL, output signals after system 4 appeared to be downsampled and
% %played in a lower tone as compared to the original. Among the three
% %L-values, output signals that are generated with L=400 are more audible
% %compared to the other signals generated from L=50 and L=100._

%% Is the system BIBO stable?
%%
% _System 1 is a BIBO stable system as poles are inside the
% boundary of the system as shown in the zplane. Also, impulse
% response graph of system 1 shows that it is finite thus, the values are
% bounded and does not change in value._

as1 = [1]; % system 1 output coef
bs1 = [0.5 0.5]; % system 1 input coef
figure(); zplane(bs1,as1); % generate z-pole of system 1;
s1N=1000;
s1n=0:s1N-1;
s1x = (s1n==0);
s1y=filter(bs1,as1,s1x);
figure(); stem(s1n,s1y);

%%
% _System 2 is NOT a BIBO stable system with its poles are outside of
% the system's circle as shown in the zplane. This is supported by
% the impulse response graph with its amplitude changes by near end of time._

as2 = [1 2 2]; % system 2 output coef
bs2 = [1]; % system 2 input coef
figure(); zplane(bs2,as2); % generate z-pole of system 2;
s2N=1000;
s2n=0:s2N-1;
s2x = (s2n==0);
s2y=filter(bs2,as2,s2x);
figure(); stem(s2n,s2y);

%%
% _System 3 is NOT a BIBO stable system as its poles are outside
% of the system's circle as shown in the zplane. Similar to system 2,
% its impulse reponse graph shows that its amplitude chanegs by near end of time._

as3 = [1 2 2]; % system 3 output coef
bs3 = [1.5 0.5]; % system 3 input coef
figure(); zplane(bs3,as3); % generate z-pole of system 3;
s3N=1000;
s3n=0:s3N-1;
s3x = (s3n==0);
s3y=filter(bs3,as3,s3x);
figure(); stem(s3n,s3y);

%%
% _System 4 is a BIBO stable system as poles are present inside
% the circle of the system as shown in the zplane below. Additionally,
% its impulse response grah shows that its amplitude remains constant until the end._

as4 = [1 -0.5 -0.5]; % system 3 output coef
bs4 = [1]; % system 3 input coef
figure();
zplane(bs4,as4); % generate z-pole of system 4;
s4N=1000;
s4n=0:s4N-1;
s4x = (s4n==0);
s4y=filter(bs4,as4,s4x);
figure(); stem(s4n,s4y);

%% Is the system causal?
%%
% _System 1 is a causal system as its ouputs depend on the current and
% previous inputs. For reference, consider the provided function which was
% implemented in matlab below. The values of y(n) depends on the terminal
% values of 1st term and 2nd term. Likewise, the 2nd term depends on the
% previous output of the 1st term in the inputs side._

% function y = dt_1(x)
%     % x must be defined in the main function
%     y = zeros(1,length(x)); %initialize output signal y
%     for n = 1:length(x) %5 time indices
%         if n<2
%             y(n) = x(n);
%         else
%             y(n) = 0.5*x(n) + 0.5*x(n-1);
%         end
%     end
% end

%%
% _System 2 is a causal system as its output depend only on the past,
% and present values of input and previous outputs. For reference, 
% consider the function which was implemented in matlab below.
% function y = dt_2(x)_ 

%     % x must be defined in the main function
%     y = zeros(1,length(x)); %initialize output signal y
%     for n = 1:length(x) %5 time indices
%          if n==1
%             y(n) = x(n);
%          elseif n==2
%             y(n) = x(n) - 2*y(n-1);
%          else
%             y(n) = x(n) - 2*y(n-1)- 2*y(n-2);
%          end
%     end
% end

%%
% _System 3 is a causal system as its output depend only on the past,
% and present inputs and past outputs. For reference, 
% consider system 3 function with n = 1.
% present output = present inpt + past inp - previous outputs_         

% y[1] = 1.5*x[1] + 0.5*x[1-1] - 2*y[1-1] - 2*y[1-2]
% y[1] = 1.58x[1] + 0.5*x[0] - 2*y[0] - 2*y[-1]

%%
% _System 4 is also a causal system with its input is dependent on current
% and past inputs and past outputs. Consider substituting L=100, n=1 to the
% function below._

% y[1] = x[1]+0.5*y[1-100]+0.5*y[1-100-1]
% y[1] = x[1]+0.5*y[-99]+0.5*y[-100]

%% Is the system FIR or IIR?
%%
% _System 1: Since system 1 is causal and its pole is located at the
% origin, system 1 is an FIR. See the Figure output of line 418 in the 
% attached file below for evidence. As a rule of thumb, if the system is
% causal and its poles are located at the origin, then it is FIR._

%%
% _System 2: Since system 2 is causal and its poles are beyond the boundary
% of the system nor the origin, it is an IIR. See the figure output of line 430 in the
% ttached file below for evidence._

%%
% _System 3: Since system 3 is causal and its poles are outside the boundary
% of the system nor the origin, it is an IIR. See the figure output of line
% 442 in the attached file below._

%%
% _System 4: Since system 4 is causal and its poles are not in the origin
% but within the boundary of the system, it is an IIR. See the figure 
% output of line 452 in the attached file below._

%% What does the system do?

% System 1
figure(); freqz(bs1,as1);
%%
% _OBSERVATION: System 1 shifted the input signal from 0 to -50dB with a
% phase of 0 to -90 degrees per cycle. This results for the output signal
% to be audible and be sense to have been no difference with the original
% audio signal._

% System 2
figure(); freqz(bs2,as2);
%%
% _OBSERVATION: Sysem 2 shifted the magnitude of the input signal from -14
% to 3dB with a phase of 0 to 350 degrees per cycle. This makes the output
% signal of the audio to be inaudible on human-level._

% System 3
figure(); freqz(bs3,as3);
%%
% _OBSERVATION: System 3 shifted the magnitude of the input signal from -8dB
% %to 5dB with a phase of 0 to 360 degrees per cycle. This results for the
% %audio output signal to be audible just like system 1 and be sensed like
% %there is no difference with the original audio due to small difference._

% System 4 
figure(); freqz(bs4,as4);
%%
% _OBSERVATION: System 4 shifted the magnitude of the input signal
% logarithmically from 40 to -3dB with a phase of -90 to 0 degrees per
% cycle. This makes the output signal to be audible on human level but with
% a mixed of a white-noise-like-audio or downsampled signal._

%% Overall
%%
% _System 1 is an FIR, a causal, and BIBO stable system that shifted the input
% audio from 0 to -50dB with a phase of 0 to -90 degrees per cycle and
% resulted to an output audio signal that is audible to human-level and
% having a little to no difference quality with the original audio signals._

%%
% _System 2 is an IIR, causal, and BIBO unstable system that shifted the
% input audio signals from -14 to 3dB with a phase of 0-350 degrees per
% cycle and makes the resulting output signal inaudible on human-level._

%%
% _System 3 is an IIR, causal, and BIBO unstable system that shifted the
% input audio signals from -8 to 5dB with a phase of 0-360 degrees per
% cycle and outcome an output signal that is audible to the human-level
% with quality as system 1 and the original audio signal._

%%
% _System 4 on the other hand is an IIR, causal, and BIBO stable system that
% log-shifted input audio signal from 40 to -3dB with a phase angle of -90 to 0
% degrees making the output signal to be audible on human-level but with
% white-noise-sound-included or of audible lower quality as the original audio
% signal._ 