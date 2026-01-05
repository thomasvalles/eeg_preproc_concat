# RF Steps Read-Me
EEGLAB PLUGINS:
- clean_rawdata
- firfilt
- ICLabel
- dipfit
- ANTeepimport
- FASTER

## Step0_SetDirectories.m

- work_dir is contains a folder for each interrogation
- subjid is the name of the folder for the interrogation you want to process
- [work_dir subjid] should contain a folder called "Step1_InputCNT", which holds the raw CNT from the interrogation

## Step1_Merge.m

- **Gets all the .cnt files from that date for that subject, merges into a single .set file**
- Note that .cnt files are ordered for merging based on time stamp automatically by matlab.
- For files with new naming convention, pay attention to the ordering to make sure the EEG are getting merged in correct order, since they are ordered automatically with the dir command

  ### CPz

  - If CPz doesn't exist, make CPz channel for ground (replace EOG label)
  - If channel 32 is not all zeros, sets all vals equal to zero

## Step2_BrushAndCut.m

- **Detects pulses so we can cut up the EEG around the pulses**
- If there are triggers, just uses triggers to find pulses, cuts around the trains.
- If there are no triggers, use one of three methods to find the pulses using peak detection. (Note-- if we do some temporary filtering here that should fix a lot of issues, don't know why we didn't do before.)
- Calculates stimulation frequency of train based on time between pulses
- After that:

  1.  cuts out artifact (.005s pre, .020s post)
  2.  Pre and Post are separately demeaned
  3.  Detrended
  4.  Interpolates around the pulse in the cutout section (.025s in total)
  5.  If the variance of 100ms of the interpolated section is too high, add some gaussian noise (note that I don't exactly get why this is here)

  --> **Outputs 1. triggers .m data file with frequency information (order of epochs), 2. an EEG with epochs of 1s pre train, 2s post train for each train (epoch)**

## Step3_FASTERauto.m

- Loads a FASTER job file with the following settings:

  1.  High pass 1Hz
  2.  Low Pass filter 55.0 Hz
  3.  Downsample to 1000Hz
  4.  Run ICA
  5.  Reject bad components
  6.  Interpolate bad channels

  --> **Faster Out EEG files**

## Step4_EntrainmentVIHalfHz.m

- Visually inspect and mark bad epochs
  Note: Epochs aren't actually removed anymore. We just mark them and then make any corresponding values as NaN in analyses. This makes it easier to keep track of.

  --> **New Trigger File with bad epochs marked**

## Step5_Entrainment_ASC_InputOutput_WholeBrain_Gabor.m

- Calculates SCC and Power via Gabor transform.
- Sets time of interest as 1s Pre and only 1s Post (note we have 2s of Post data)
- For each epoch (stimulation train):
  - calculate power using Gabor transform
  - given absolute power, for each channel, calculate SCC with F3 as seed channel.
    - Do this at ever possible frequency from the Gabor transform
    - eg calculate SCC from F3 <-> Fp1 , with "center frequency" 2 Hz, 2.1 Hz, 2.2 Hz.... 20 Hz etc
- SCC Calculation between F3 and chan_2 with $f$ center frequency:

$$ SCC(f, F3,chan_2) = corr(PSD_{F3}^{[f-2,f+2]},PSD_{chan_2}^{[f-2,f+2]})$$


Where $PSD_{F3}, PSD_{chan_2}$ are **power spectra** of respective channels. Note that SCC calculation at frequency $f$ looks at the correlation between two power spectra in a 4Hz neighborhood($\pm$ 2 Hz) of $f$.

- For any epochs that were marked in visual inspection (step 6), make the SCC or Power values all NaN for this train (epoch).

  --> **Save data files in various formats**

## Step6_FrequencyDecision_UpdatedModel_Repeats.m

- This for the matched SCC model (2022-ish). There's also the 13-17 Hz model written in as another option.

- Channels in ROI : Fpz, AF4 (seed: F3)
  - we will eventually take **median** over the ROI centered at these two channels, according to closest channels on EEG headmap.
- For each unique frequency (eg $f=10.5Hz$ stimulation), consider "matched" window, ie 8.5-12.5Hz. Suppose there are 4 epochs of $f$. Then can calculate the average change in SCC in the $[8.5,12.5]$ band for each epoch (ie 4x64 array), considering all possible channels.
- Then, take the median value over each ROI (ie median of the change at Fpz and each of its neighbors), so that there is 1 value per ROI per trial.
- For the boxplots ,we can stop simply multiply each of these trials by the regression weights.
- For a single "median" value for a stimulation frequency, take the median value over the 4 epochs (eg 1 value per epoch per ROI --> 1 value per ROI). (This is how the "training" for the model was done)
- Then multiply the value by the regression weights to get a "connectivity change" value for that particular frequency.

**Important notes** :
In this script, "matched SCC" is not just a single SCC value calculated in $[stim freq -2 , stim freq + 2]$ band (ie centered at stim freq). It is an **average** of SCC values in this band -- ie: $$mean\{SCC_{stim freq -2}, SCC_{stim freq -1.9}, ..., SCC_{stim freq},..., SCC_{stim freq + 1.9}, SCC_{stim freq + 2}\}$$
Rather than simply

$$ SCC_{stim freq}$$

In other words, an individual SCC value is calculated as the PSD correlation in a $\pm 2$ Hz band, **then** second, we average many SCC values over a $\pm 2$ Hz band.
