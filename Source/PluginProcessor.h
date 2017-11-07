/*
  This code accompanies the textbook:
 
  Digital Audio Effects: Theory, Implementation and Application
  Joshua D. Reiss and Andrew P. McPherson
 
  ---
 
  Compressor: dynamic range compression effect
  See textbook Chapter 6: Dynamics Processing
 
  Code by Joshua Reiss, Brecht de Man and Andrew McPherson
 
  ---

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __PLUGINPROCESSOR_H_88534BAA__
#define __PLUGINPROCESSOR_H_88534BAA__

#include "../JuceLibraryCode/JuceHeader.h"
#include <math.h>
class CompressorAudioProcessor  : public AudioProcessor
{
public:
    CompressorAudioProcessor();
    ~CompressorAudioProcessor();
	
	int bufferSize;
	
    void prepareToPlay (double sampleRate, int samplesPerBlock);
    void releaseResources();
	void processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages);
	
	void compressorLF(AudioSampleBuffer &buffer, int m);// compressor functions
    void compressorMF(AudioSampleBuffer &buffer, int m);// compressor functions
    void compressorHF(AudioSampleBuffer &buffer, int m);// compressor functions

	template <class T> const T& max ( const T& a, const T& b );

    AudioProcessorEditor* createEditor();

    bool hasEditor() const;

	AudioPlayHead::CurrentPositionInfo lastPosInfo;
 
	int round(float inn);
    const String getName() const;

    int getNumParameters();

    float getParameter (int index);
    void setParameter (int index, float newValue);

    const String getParameterName (int index);
    const String getParameterText (int index);

    const String getInputChannelName (int channelIndex) const;
    const String getOutputChannelName (int channelIndex) const;
    bool isInputChannelStereoPair (int index) const;
    bool isOutputChannelStereoPair (int index) const;

    bool silenceInProducesSilenceOut() const;
    double getTailLengthSeconds() const;
    bool acceptsMidi() const;
    bool producesMidi() const;

    int getNumPrograms();
    int getCurrentProgram();
    void setCurrentProgram (int index);
    const String getProgramName (int index);
    void changeProgramName (int index, const String& newName);

    //==============================================================================
    void getStateInformation (MemoryBlock& destData);
    void setStateInformation (const void* data, int sizeInBytes);

	float getThresholdLF();
	float getRatioLF();
	float getGainLF();
	float getAttackTimeLF();
	float getReleaseTimeLF();
	void setThresholdLF(float Tlf);
	void setGainLF(float Glf);
	void setRatioLF(float Rlf);
	void setAttackTimeLF(float Alf);
	void setReleaseTimeLF(float Rlf);
    
    float getThresholdMF();
    float getRatioMF();
    float getGainMF();
    float getAttackTimeMF();
    float getReleaseTimeMF();
    void setThresholdMF(float Tmf);
    void setGainMF(float Gmf);
    void setRatioMF(float Rmf);
    void setAttackTimeMF(float Amf);
    void setReleaseTimeMF(float Rmf);
    
    float getThresholdHF();
    float getRatioHF();
    float getGainHF();
    float getAttackTimeHF();
    float getReleaseTimeHF();
    void setThresholdHF(float Thf);
    void setGainHF(float Ghf);
    void setRatioHF(float Rhf);
    void setAttackTimeHF(float Ahf);
    void setReleaseTimeHF(float Rhf);
    
	void resetAll();
    

	// parameters

	bool compressorONOFF;
	int M;
	bool autoTime;
    
    //FILTER STUFF
    IIRFilter lowPassL1, lowPassL2, lowPassR1, lowPassR2; //Filters for low band
    IIRFilter lowBandPassL1, lowBandPassL2, lowBandPassR1, lowBandPassR2; //Filters for low band pass
    IIRFilter highBandPassL1, highBandPassL2, highBandPassR1, highBandPassR2; //Filters for high band pass
    IIRFilter highPassL1, highPassL2, highPassR1, highPassR2; //Filters for high pass
    
    IIRCoefficients coeff; //Coefficient variable
    
    
    //Filter Parameters
    float getLCF();
    float getHCF();
    void setLCF(float fcLow);
    void setHCF(float fcHigh);
    
    //Cuttoff frequency for change detection
    float lastfcLow;
    float lastfcHigh;
    
    //Initial cuttoff frequencies
    float fcLow;
    float fcHigh;
    
    float LCF;
    float HCF;

private:
    AudioSampleBuffer inputBuffer;

//	int bufferSize;
    //these are used to persist UI's size- values are stored along with filter's other parameters, and UI component will update them when it gets resized.
	int lastUIWidth, lastUIHeight;
		
	HeapBlock <float> x_g_lf, x_l_lf,y_g_lf, y_l_lf,c_lf;// input, output, control
    HeapBlock <float> x_g_mf, x_l_mf,y_g_mf, y_l_mf,c_mf;// input, output, control
    HeapBlock <float> x_g_hf, x_l_hf,y_g_hf, y_l_hf,c_hf;// input, output, control
    
    
		// Compressor parameters for each frequency band
	float ratioLF,thresholdLF,makeUpGainLF,tauAttackLF,tauReleaseLF,alphaAttackLF,alphaReleaseLF,yL_prevLF;
    float ratioHF,thresholdHF,makeUpGainHF,tauAttackHF,tauReleaseHF,alphaAttackHF,alphaReleaseHF,yL_prevHF;
    float ratioMF,thresholdMF,makeUpGainMF,tauAttackMF,tauReleaseMF,alphaAttackMF,alphaReleaseMF,yL_prevMF;


	int nhost;
	int samplerate;

	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (CompressorAudioProcessor);
};

#endif  // __PLUGINPROCESSOR_H_88534BAA__
