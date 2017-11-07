/*
  ==============================================================================

  This is an automatically generated GUI class created by the Introjucer!

  Be careful when adding custom code to these files, as only the code within
  the "//[xyz]" and "//[/xyz]" sections will be retained when the file is loaded
  and re-saved.

  Created with Introjucer version: 4.1.0

  ------------------------------------------------------------------------------

  The Introjucer is part of the JUCE library - "Jules' Utility Class Extensions"
  Copyright (c) 2015 - ROLI Ltd.

  ==============================================================================
*/

#ifndef __JUCE_HEADER_388D2EBDBB646058__
#define __JUCE_HEADER_388D2EBDBB646058__

//[Headers]     -- You can add your own extra header files here --
//[/Headers]



//==============================================================================
/**
                                                                    //[Comments]
    An auto-generated component, created by the Jucer.

    Describe your class and how it works here!
                                                                    //[/Comments]
*/
class CompressorAudioProcessorEditor  : public AudioProcessorEditor,
                                        public Timer,
                                        public ButtonListener,
                                        public SliderListener
{
public:
    //==============================================================================
    CompressorAudioProcessorEditor (CompressorAudioProcessor* ownerFilter);
    ~CompressorAudioProcessorEditor();

    //==============================================================================
    //[UserMethods]     -- You can add your own custom methods in this section.
	void timerCallback();
    //[/UserMethods]

    void paint (Graphics& g);
    void resized();
    void buttonClicked (Button* buttonThatWasClicked);
    void sliderValueChanged (Slider* sliderThatWasMoved);

    // Binary resources:
    static const char* brushedMetalDark_jpg;
    static const int brushedMetalDark_jpgSize;
    static const char* c4dm_png2;
    static const int c4dm_png2Size;
    static const char* qmul_png2;
    static const int qmul_png2Size;
    static const char* knobstrip_png;
    static const int knobstrip_pngSize;
    static const char* scaleLr_png;
    static const int scaleLr_pngSize;


private:
    //[UserVariables]   -- You can add your own custom variables in this section.

    ScopedPointer<ResizableCornerComponent> resizer;
    ComponentBoundsConstrainer resizeLimits;



	AudioPlayHead::CurrentPositionInfo lastDisplayedPosition;

    CompressorAudioProcessor* getProcessor() const
    {
        return static_cast <CompressorAudioProcessor*> (getAudioProcessor());
    }

    void displayPositionInfo (const AudioPlayHead::CurrentPositionInfo& pos);

    //[/UserVariables]

    //==============================================================================
    ScopedPointer<TextButton> buttonONOFF;
    ScopedPointer<Label> title;
    ScopedPointer<Slider> sliderThresholdLF;
    ScopedPointer<Label> labelThresholdLF;
    ScopedPointer<Slider> sliderRatioLF;
    ScopedPointer<Label> labelRatioLF;
    ScopedPointer<Slider> sliderGainLF;
    ScopedPointer<Label> labelGainLF;
    ScopedPointer<Slider> sliderAttackLF;
    ScopedPointer<Label> labelAttackLF;
    ScopedPointer<Slider> sliderReleaseLF;
    ScopedPointer<Label> labelReleaseLF;
    ScopedPointer<Slider> sliderThresholdMF;
    ScopedPointer<Label> labelThresholdMF;
    ScopedPointer<Slider> sliderRatioMF;
    ScopedPointer<Label> labelRatioMF;
    ScopedPointer<Slider> sliderGainMF;
    ScopedPointer<Label> labelGainMF;
    ScopedPointer<Slider> sliderAttackMF;
    ScopedPointer<Label> labelAttackMF;
    ScopedPointer<Slider> sliderReleaseMF;
    ScopedPointer<Label> labelReleaseMF;
    ScopedPointer<Slider> sliderThresholdHF;
    ScopedPointer<Label> labelThresholdHF;
    ScopedPointer<Slider> sliderRatioHF;
    ScopedPointer<Label> labelRatioHF;
    ScopedPointer<Slider> sliderGainHF;
    ScopedPointer<Label> labelGainHF;
    ScopedPointer<Slider> sliderAttackHF;
    ScopedPointer<Label> labelAttackHF;
    ScopedPointer<Slider> sliderReleaseHF;
    ScopedPointer<Label> labelReleaseHF;
    ScopedPointer<Slider> LCF;
    ScopedPointer<Slider> HCF;


    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (CompressorAudioProcessorEditor)
};

//[EndFile] You can add extra defines here...
//[/EndFile]

#endif   // __JUCE_HEADER_388D2EBDBB646058__
