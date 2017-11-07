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

#include "PluginProcessor.h"
#include "PluginEditor.h"
CompressorAudioProcessor::CompressorAudioProcessor()
	// Initializer List
	:
	inputBuffer(1,1),
	nhost(0)
{
	lastUIWidth = 850; 
    lastUIHeight = 650;
    lastPosInfo.resetToDefault();
}
CompressorAudioProcessor::~CompressorAudioProcessor()
{
}
//==============================================================================
void CompressorAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback initialisation that you need.
	M = round(getTotalNumInputChannels()/2);
	samplerate = (float)getSampleRate();
	bufferSize = getBlockSize();
	// Allocate a lot of dynamic memory here for low band
	x_g_lf					.allocate(bufferSize, true);
	x_l_lf					.allocate(bufferSize, true);
	y_g_lf					.allocate(bufferSize, true);
	y_l_lf					.allocate(bufferSize, true);
	c_lf					.allocate(bufferSize, true);
    
    // Allocate a lot of dynamic memory here for mid band
    x_g_mf					.allocate(bufferSize, true);
    x_l_mf					.allocate(bufferSize, true);
    y_g_mf					.allocate(bufferSize, true);
    y_l_mf					.allocate(bufferSize, true);
    c_mf					.allocate(bufferSize, true);
    
    // Allocate a lot of dynamic memory here for high band
    x_g_hf					.allocate(bufferSize, true);
    x_l_hf					.allocate(bufferSize, true);
    y_g_hf					.allocate(bufferSize, true);
    y_l_hf					.allocate(bufferSize, true);
    c_hf					.allocate(bufferSize, true);
    
    
	yL_prevLF=0;
    yL_prevHF=0;
    yL_prevMF=0;

    
	autoTime = false;
	compressorONOFF = false;
	resetAll();
    
    
     fcLow=500;
     fcHigh=5000;
     lastfcLow=500;
     lastfcHigh=5000;
    
    //===================================================//
    //          Calculate Coefficients                   //
    //===================================================//
    
    lowPassL1.setCoefficients(coeff.makeLowPass(sampleRate, fcLow));
    lowPassL2.setCoefficients(coeff.makeLowPass(sampleRate, fcLow));
    lowPassR1.setCoefficients(coeff.makeLowPass(sampleRate, fcLow));
    lowPassR2.setCoefficients(coeff.makeLowPass(sampleRate, fcLow));
    
    lowBandPassL1.setCoefficients(coeff.makeHighPass(sampleRate, fcLow));
    lowBandPassL2.setCoefficients(coeff.makeHighPass(sampleRate, fcLow));
    lowBandPassR1.setCoefficients(coeff.makeHighPass(sampleRate, fcLow));
    lowBandPassR2.setCoefficients(coeff.makeHighPass(sampleRate, fcLow));
    
    highBandPassL1.setCoefficients(coeff.makeLowPass(sampleRate, fcHigh));
    highBandPassL2.setCoefficients(coeff.makeLowPass(sampleRate, fcHigh));
    highBandPassR1.setCoefficients(coeff.makeLowPass(sampleRate, fcHigh));
    highBandPassR2.setCoefficients(coeff.makeLowPass(sampleRate, fcHigh));
    
    highPassL1.setCoefficients(coeff.makeHighPass(sampleRate, fcHigh));
    highPassL2.setCoefficients(coeff.makeHighPass(sampleRate, fcHigh));
    highPassR1.setCoefficients(coeff.makeHighPass(sampleRate, fcHigh));
    highPassR2.setCoefficients(coeff.makeHighPass(sampleRate, fcHigh));
    
    
    
    
}
void CompressorAudioProcessor::releaseResources()
{
    // When playback stops, you can use this to free up any spare memory, etc.
}
void CompressorAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
    
    int numChannels=buffer.getNumChannels();
    int numSamples=buffer.getNumSamples();
    
    
    //===================================================//
    //   Re-Calculate Coefficients if Cuttoff Changes    //
    //===================================================//
    
    if (fcLow != lastfcLow) {
        lowPassL1.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        lowPassL2.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        lowPassR1.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        lowPassR2.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        
        lowBandPassL1.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lowBandPassL2.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lowBandPassR1.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lowBandPassR2.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lastfcLow=fcLow;
    }
    
    if (fcHigh != lastfcHigh) {
        lowPassL1.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        lowPassL2.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        lowPassR1.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        lowPassR2.setCoefficients(coeff.makeLowPass(samplerate, fcLow));
        
        lowBandPassL1.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lowBandPassL2.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lowBandPassR1.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lowBandPassR2.setCoefficients(coeff.makeHighPass(samplerate, fcLow));
        lastfcHigh=fcHigh;
    }
    
    
    //===================================================//
    //   Create and Empty Buffers                        //
    //===================================================//
    
    AudioSampleBuffer lpFilterBuffer;
    AudioSampleBuffer hpFilterBuffer;
    AudioSampleBuffer bpFilterBuffer;
    
    AudioSampleBuffer lpCompFilterBuffer;
    AudioSampleBuffer hpCompFilterBuffer;
    AudioSampleBuffer bpCompFilterBuffer;
    
    lpFilterBuffer.setSize(2, numSamples);
    lpFilterBuffer.clear();
    
    hpFilterBuffer.setSize(2, numSamples);
    hpFilterBuffer.clear();
    
    bpFilterBuffer.setSize(2, numSamples);
    bpFilterBuffer.clear();
    
    lpCompFilterBuffer.setSize(2, numSamples);
    lpCompFilterBuffer.clear();
    
    hpCompFilterBuffer.setSize(2, numSamples);
    hpCompFilterBuffer.clear();
    
    bpCompFilterBuffer.setSize(2, numSamples);
    bpCompFilterBuffer.clear();
    
    lpFilterBuffer.makeCopyOf(buffer);
    hpFilterBuffer.makeCopyOf(buffer);
    bpFilterBuffer.makeCopyOf(buffer);
    
    //===================================================//
    //   Apply Filters to Buffers                        //
    //===================================================//
    
    
    lowPassL1.processSamples(lpFilterBuffer.getWritePointer(0), bufferSize);
    lowPassR1.processSamples(lpFilterBuffer.getWritePointer(1), bufferSize);
    lowPassL2.processSamples(lpFilterBuffer.getWritePointer(0), bufferSize);
    lowPassR2.processSamples(lpFilterBuffer.getWritePointer(1), bufferSize);
    
    lowBandPassL1.processSamples(bpFilterBuffer.getWritePointer(0), bufferSize);
    lowBandPassR1.processSamples(bpFilterBuffer.getWritePointer(1), bufferSize);
    lowBandPassL2.processSamples(bpFilterBuffer.getWritePointer(0), bufferSize);
    lowBandPassR2.processSamples(bpFilterBuffer.getWritePointer(1), bufferSize);
    
    
    highBandPassL1.processSamples(bpFilterBuffer.getWritePointer(0), bufferSize);
    highBandPassR1.processSamples(bpFilterBuffer.getWritePointer(1), bufferSize);
    highBandPassL2.processSamples(bpFilterBuffer.getWritePointer(0), bufferSize);
    highBandPassR2.processSamples(bpFilterBuffer.getWritePointer(1), bufferSize);
    
    
    highPassL1.processSamples(hpFilterBuffer.getWritePointer(0), bufferSize);
    highPassR1.processSamples(hpFilterBuffer.getWritePointer(1), bufferSize);
    highPassL2.processSamples(hpFilterBuffer.getWritePointer(0), bufferSize);
    highPassR2.processSamples(hpFilterBuffer.getWritePointer(1), bufferSize);
    
      
    //===================================================//
    //   Apply Compressors to Buffers                    //
    //===================================================//
    
        if (compressorONOFF)
        {
            
            for (int m = 0 ; m < M ; ++m)
            {
                if ( (thresholdLF< 0) )
                {
                    lpCompFilterBuffer.addFrom(m,0,lpFilterBuffer,m*2,0,bufferSize,0.5);
                    lpCompFilterBuffer.addFrom(m,0,lpFilterBuffer,m*2+1,0,bufferSize,0.5);
                    //compression : calculates the control voltage
                    compressorLF(inputBuffer,m); //ADD BAND NUMBER INSTEAD OF 1
                    // apply control voltage to the audio signal
                    for (int i = 0 ; i < bufferSize ; ++i)
                    {
                        lpFilterBuffer.getWritePointer(2*m+0)[i] *= c_lf[i];
                        lpFilterBuffer.getWritePointer(2*m+1)[i] *= c_lf[i];
                    }
                    lpCompFilterBuffer.clear(m,0,bufferSize);
                }
                
                if ( (thresholdMF< 0) )
                {
                    bpCompFilterBuffer.addFrom(m,0,bpFilterBuffer,m*2,0,bufferSize,0.5);
                    bpCompFilterBuffer.addFrom(m,0,bpFilterBuffer,m*2+1,0,bufferSize,0.5);
                    //compression : calculates the control voltage
                    compressorMF(bpCompFilterBuffer,m); //ADD BAND
                    // apply control voltage to the audio signal
                    for (int i = 0 ; i < bufferSize ; ++i)
                    {
                        bpFilterBuffer.getWritePointer(2*m+0)[i] *= c_mf[i];
                        bpFilterBuffer.getWritePointer(2*m+1)[i] *= c_mf[i];
                    }
                    bpCompFilterBuffer.clear(m,0,bufferSize);
                }
                if ( (thresholdHF< 0) )
                {
                    hpCompFilterBuffer.addFrom(m,0,hpFilterBuffer,m*2,0,bufferSize,0.5);
                    hpCompFilterBuffer.addFrom(m,0,hpFilterBuffer,m*2+1,0,bufferSize,0.5);
                    //compression : calculates the control voltage
                    compressorHF(hpCompFilterBuffer,m); //ADD BAND NUMBER
                    // apply control voltage to the audio signal
                    for (int i = 0 ; i < bufferSize ; ++i)
                    {
                        hpFilterBuffer.getWritePointer(2*m+0)[i] *= c_hf[i];
                        hpFilterBuffer.getWritePointer(2*m+1)[i] *= c_hf[i];
                    }
                    hpCompFilterBuffer.clear(m,0,bufferSize);
                }
                
            }
            for(int channel=0; channel<numChannels; channel++){
                for(int samples=0; samples<numSamples; samples++){
                    buffer.copyFrom(channel, 0, lpFilterBuffer, channel, 0, numSamples);
                    buffer.addFrom(channel, 0, bpFilterBuffer, channel, 0, numSamples);
                    buffer.addFrom(channel, 0, hpFilterBuffer, channel, 0, numSamples);
                }
            }
        }
	}

// Low compressor function
void CompressorAudioProcessor::compressorLF(AudioSampleBuffer &buffer, int m)
{
	alphaAttackLF = exp(-1/(0.001 * samplerate * tauAttackLF));
	alphaReleaseLF= exp(-1/(0.001 * samplerate * tauReleaseLF));
	for (int i = 0 ; i < bufferSize ; ++i)
	{
		//Level detection- estimate level using peak detector
		if (fabs(buffer.getWritePointer(m)[i]) < 0.000001) x_g_lf[i] =-120;
		else x_g_lf[i] =20*log10(fabs(buffer.getWritePointer(m)[i]));
		//Gain computer- static apply input/output curve
		if (x_g_lf[i] >= thresholdLF) y_g_lf[i] = thresholdLF+ (x_g_lf[i] - thresholdLF) / ratioLF;
		else y_g_lf[i] = x_g_lf[i];
		x_l_lf[i] = x_g_lf[i] - y_g_lf[i];
		//Ballistics- smoothing of the gain 
		if (x_l_lf[i]>yL_prevLF)  y_l_lf[i]=alphaAttackLF * yL_prevLF+(1 - alphaAttackLF) * x_l_lf[i] ;
		else				 y_l_lf[i]=alphaReleaseLF* yL_prevLF+(1 - alphaReleaseLF) * x_l_lf[i] ;
		//find control
		c_lf[i] = pow(10,(makeUpGainLF - y_l_lf[i])/20);
		yL_prevLF=y_l_lf[i];
	}
}

// Mid compressor function
void CompressorAudioProcessor::compressorMF(AudioSampleBuffer &buffer, int m)
{
    alphaAttackMF = exp(-1/(0.001 * samplerate * tauAttackMF));
    alphaReleaseMF= exp(-1/(0.001 * samplerate * tauReleaseMF));
    for (int i = 0 ; i < bufferSize ; ++i)
    {
        //Level detection- estimate level using peak detector
        if (fabs(buffer.getWritePointer(m)[i]) < 0.000001) x_g_mf[i] =-120;
        else x_g_mf[i] =20*log10(fabs(buffer.getWritePointer(m)[i]));
        //Gain computer- static apply input/output curve
        if (x_g_mf[i] >= thresholdMF) y_g_mf[i] = thresholdMF+ (x_g_mf[i] - thresholdMF) / ratioMF;
        else y_g_mf[i] = x_g_mf[i];
        x_l_mf[i] = x_g_mf[i] - y_g_mf[i];
        //Ballistics- smoothing of the gain
        if (x_l_mf[i]>yL_prevMF)  y_l_mf[i]=alphaAttackMF * yL_prevMF+(1 - alphaAttackMF) * x_l_mf[i] ;
        else				 y_l_mf[i]=alphaReleaseMF* yL_prevMF+(1 - alphaReleaseMF) * x_l_mf[i] ;
        //find control
        c_mf[i] = pow(10,(makeUpGainMF - y_l_mf[i])/20);
        yL_prevMF=y_l_mf[i];
    }
}

//High compressor function
void CompressorAudioProcessor::compressorHF(AudioSampleBuffer &buffer, int m)
{
    alphaAttackHF = exp(-1/(0.001 * samplerate * tauAttackHF));
    alphaReleaseHF= exp(-1/(0.001 * samplerate * tauReleaseHF));
    for (int i = 0 ; i < bufferSize ; ++i)
    {
        //Level detection- estimate level using peak detector
        if (fabs(buffer.getWritePointer(m)[i]) < 0.000001) x_g_hf[i] =-120;
        else x_g_hf[i] =20*log10(fabs(buffer.getWritePointer(m)[i]));
        //Gain computer- static apply input/output curve
        if (x_g_hf[i] >= thresholdHF) y_g_hf[i] = thresholdHF+ (x_g_hf[i] - thresholdHF) / ratioHF;
        else y_g_hf[i] = x_g_hf[i];
        x_l_hf[i] = x_g_hf[i] - y_g_hf[i];
        //Ballistics- smoothing of the gain
        if (x_l_hf[i]>yL_prevHF)  y_l_hf[i]=alphaAttackHF * yL_prevHF+(1 - alphaAttackHF ) * x_l_hf[i] ;
        else				 y_l_hf[i]=alphaReleaseHF* yL_prevHF+(1 - alphaReleaseHF) * x_l_hf[i] ;
        //find control
        c_hf[i] = pow(10,(makeUpGainHF - y_l_hf[i])/20);
        yL_prevHF=y_l_hf[i];
    }
}
template <class T> const T& CompressorAudioProcessor::max( const T& a, const T& b )
{
  return (a < b) ? b : a;
}
void CompressorAudioProcessor::resetAll()
{
    tauAttackLF=0;tauReleaseLF = 0;
    alphaAttackLF=0;alphaReleaseLF = 0;
    thresholdLF = 0;
    ratioLF= 1;
    makeUpGainLF= 0;
    yL_prevLF=0;
    
    tauAttackMF=0;tauReleaseMF = 0;
    alphaAttackMF=0;alphaReleaseMF = 0;
    thresholdMF = 0;
    ratioMF= 1;
    makeUpGainMF= 0;
    yL_prevMF=0;
    
    tauAttackHF=0;tauReleaseHF = 0;
    alphaAttackHF=0;alphaReleaseHF = 0;
    thresholdHF = 0;
    ratioHF= 1;
    makeUpGainHF= 0;
    yL_prevHF=0;
    
	for (int i = 0 ; i < bufferSize ; ++i)
	{
		x_g_lf[i] = 0;	y_g_lf[i] = 0;
		x_l_lf[i] = 0;	y_l_lf[i] = 0;
		c_lf[i] = 0;
	}
    for (int i = 0 ; i < bufferSize ; ++i)
    {
        x_g_mf[i] = 0;	y_g_mf[i] = 0;
        x_l_mf[i] = 0;	y_l_mf[i] = 0;
        c_mf[i] = 0;
    }
    for (int i = 0 ; i < bufferSize ; ++i)
    {
        x_g_hf[i] = 0;	y_g_hf[i] = 0;
        x_l_hf[i] = 0;	y_l_hf[i] = 0;
        c_hf[i] = 0;
    }
}

//===================================================//
//          Low Comp Get & Set Methods               //
//===================================================//
//////////////////////////////////////////////
float CompressorAudioProcessor::getThresholdLF()
{
	return thresholdLF;
}
float CompressorAudioProcessor::getRatioLF()
{
	return ratioLF;
}
float CompressorAudioProcessor::getGainLF()
{
	return makeUpGainLF;//problem?
}
float CompressorAudioProcessor::getAttackTimeLF()
{
	return tauAttackLF;
}
float CompressorAudioProcessor::getReleaseTimeLF()
{
	return tauReleaseLF;
}
////////////////////////////////////////////////////////
void CompressorAudioProcessor::setThresholdLF(float Tlf)
{
	thresholdLF= Tlf;
}
void CompressorAudioProcessor::setGainLF(float Glf)
{
	makeUpGainLF= Glf;
}
void CompressorAudioProcessor::setRatioLF(float Rlf)
{
	ratioLF= Rlf;
}
void CompressorAudioProcessor::setAttackTimeLF(float Alf)
{
	tauAttackLF = Alf;
}
void CompressorAudioProcessor::setReleaseTimeLF(float Rlf)
{
	tauReleaseLF = Rlf;
}

//===================================================//
//          MID Comp Get & Set Methods               //
//===================================================//
//////////////////////////////////////////////
float CompressorAudioProcessor::getThresholdMF()
{
    return thresholdMF;
}
float CompressorAudioProcessor::getRatioMF()
{
    return ratioMF;
}
float CompressorAudioProcessor::getGainMF()
{
    return makeUpGainMF;//problem?
}
float CompressorAudioProcessor::getAttackTimeMF()
{
    return tauAttackMF;
}
float CompressorAudioProcessor::getReleaseTimeMF()
{
    return tauReleaseMF;
}
////////////////////////////////////////////////////////
void CompressorAudioProcessor::setThresholdMF(float Tmf)
{
    thresholdMF= Tmf;
}
void CompressorAudioProcessor::setGainMF(float Gmf)
{
    makeUpGainMF= Gmf;
}
void CompressorAudioProcessor::setRatioMF(float Rmf)
{
    ratioMF= Rmf;
}
void CompressorAudioProcessor::setAttackTimeMF(float Amf)
{
    tauAttackMF = Amf;
}
void CompressorAudioProcessor::setReleaseTimeMF(float Rmf)
{
    tauReleaseMF = Rmf;
}

//===================================================//
//          HIGH Comp Get & Set Methods              //
//===================================================//
//////////////////////////////////////////////
float CompressorAudioProcessor::getThresholdHF()
{
    return thresholdHF;
}
float CompressorAudioProcessor::getRatioHF()
{
    return ratioHF;
}
float CompressorAudioProcessor::getGainHF()
{
    return makeUpGainHF;//problem?
}
float CompressorAudioProcessor::getAttackTimeHF()
{
    return tauAttackHF;
}
float CompressorAudioProcessor::getReleaseTimeHF()
{
    return tauReleaseHF;
}
////////////////////////////////////////////////////////
void CompressorAudioProcessor::setThresholdHF(float Thf)
{
    thresholdHF= Thf;
}
void CompressorAudioProcessor::setGainHF(float Ghf)
{
    makeUpGainHF= Ghf;
}
void CompressorAudioProcessor::setRatioHF(float Rhf)
{
    ratioHF= Rhf;
}
void CompressorAudioProcessor::setAttackTimeHF(float Ahf)
{
    tauAttackHF = Ahf;
}
void CompressorAudioProcessor::setReleaseTimeHF(float Rhf)
{
    tauReleaseHF = Rhf;
}



//===================================================//
//          Cutoff Filter Set Methods                //
//===================================================//
void CompressorAudioProcessor::setLCF(float LCF)
{
    fcLow = LCF;
}

void CompressorAudioProcessor::setHCF(float HCF)
{
    fcHigh = HCF;
}


//===================================================//
//          Cutoff Filter Get Methods                //
//===================================================//
float CompressorAudioProcessor::getLCF()
{
    return LCF;
}

float CompressorAudioProcessor::getHCF()
{
    return HCF;
    
}


//===================================================//

bool CompressorAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}
AudioProcessorEditor* CompressorAudioProcessor::createEditor()
{
    return new CompressorAudioProcessorEditor (this);
}
//==============================================================================
void CompressorAudioProcessor::getStateInformation (MemoryBlock& destData)
{
//Use this to store your parameters in memory block, either as raw data, or use XML or ValueTree classes as intermediaries to make it easy to save and load complex data.
}
void CompressorAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
// Use this to restore your parameters from this memory block, whose contents will have been created by the getStateInformation() call.
}
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new CompressorAudioProcessor();
}
int CompressorAudioProcessor::round(float inn)
{
	if (inn > 0) return (int) (inn + 0.5);
	else return (int) (inn - 0.5);
}
const String CompressorAudioProcessor::getName() const
{
    return JucePlugin_Name;
}
int CompressorAudioProcessor::getNumParameters()
{
    return 0;
}
float CompressorAudioProcessor::getParameter (int index)
{
    return 0.0f;
}
void CompressorAudioProcessor::setParameter (int index, float newValue)
{
}
const String CompressorAudioProcessor::getParameterName (int index)
{
    return String::empty;
}
const String CompressorAudioProcessor::getParameterText (int index)
{
    return String::empty;
}
const String CompressorAudioProcessor::getInputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}
const String CompressorAudioProcessor::getOutputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}
bool CompressorAudioProcessor::isInputChannelStereoPair (int index) const
{
    return true;
}
bool CompressorAudioProcessor::isOutputChannelStereoPair (int index) const
{
    return true;
}
bool CompressorAudioProcessor::silenceInProducesSilenceOut() const
{
#if JucePlugin_SilenceInProducesSilenceOut
    return true;
#else
    return false;
#endif
}

double CompressorAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}
bool CompressorAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
}
bool CompressorAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
}
int CompressorAudioProcessor::getNumPrograms()
{
    return 1;
}
int CompressorAudioProcessor::getCurrentProgram()
{
    return 0;
}
void CompressorAudioProcessor::setCurrentProgram (int index)
{
}
const String CompressorAudioProcessor::getProgramName (int index)
{
    return String::empty;
}
void CompressorAudioProcessor::changeProgramName (int index, const String& newName)
{
}