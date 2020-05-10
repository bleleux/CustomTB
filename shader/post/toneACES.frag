// Change to toneACES for the approximation that TB3 comes with
// or toneACESRef for the reference version from AMPAS' Github: https://github.com/ampas/aces-dev

#define	ToneMap	toneACESRef


vec3	toneACES( vec3 c )
{
	//ACES RRT/ODT curve fit courtesy of Stephen Hill
	vec3 a = c * (c + 0.0245786) - 0.000090537;
	vec3 b = c * (0.983729 * c + 0.4329510) + 0.238081;
	return a / b;
}

/* Begin ACES Reference code 
=============================================================================
  ACES: Academy Color Encoding System
  https://github.com/ampas/aces-dev/tree/v1.0
  
  License Terms for Academy Color Encoding System Components

  Academy Color Encoding System (ACES) software and tools are provided by the Academy under 
  the following terms and conditions: A worldwide, royalty-free, non-exclusive right to copy, modify, create
  derivatives, and use, in source and binary forms, is hereby granted, subject to acceptance of this license.

  Copyright © 2013 Academy of Motion Picture Arts and Sciences (A.M.P.A.S.). Portions contributed by
  others as indicated. All rights reserved.

  Performance of any of the aforementioned acts indicates acceptance to be bound by the following 
  terms and conditions:

   *  Copies of source code, in whole or in part, must retain the above copyright 
    notice, this list of conditions and the Disclaimer of Warranty.
   *  Use in binary form must retain the above copyright notice, this list of 
    conditions and the Disclaimer of Warranty in the documentation and/or other 
    materials provided with the distribution.
   *  Nothing in this license shall be deemed to grant any rights to trademarks,
    copyrights, patents, trade secrets or any other intellectual property of 
    A.M.P.A.S. or any contributors, except as expressly stated herein.
   *  Neither the name "A.M.P.A.S." nor the name of any other contributors to this 
    software may be used to endorse or promote products derivative of or based on
    this software without express prior written permission of A.M.P.A.S. or the
    contributors, as appropriate.

  This license shall be construed pursuant to the laws of the State of California,
  and any disputes related thereto shall be subject to the jurisdiction of the courts therein.

  Disclaimer of Warranty: THIS SOFTWARE IS PROVIDED BY A.M.P.A.S. AND CONTRIBUTORS "AS
  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND 
  NON-INFRINGEMENT ARE DISCLAIMED. IN NO EVENT SHALL A.M.P.A.S., OR ANY 
  CONTRIBUTORS OR DISTRIBUTORS, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, RESITUTIONARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  WITHOUT LIMITING THE GENERALITY OF THE FOREGOING, THE ACADEMY SPECIFICALLY 
  DISCLAIMS ANY REPRESENTATIONS OR WARRANTIES WHATSOEVER RELATED TO PATENT OR 
  OTHER INTELLECTUAL PROPERTY RIGHTS IN THE ACADEMY COLOR ENCODING SYSTEM, OR 
  APPLICATIONS THEREOF, HELD BY PARTIES OTHER THAN A.M.P.A.S.,WHETHER DISCLOSED
  OR UNDISCLOSED.
=============================================================================*/

float PI = 3.14159265359;

static const float3x3 D60_2_D65_CAT =
{
	 0.987224,   -0.00611327, 0.0159533,
	-0.00759836,  1.00186,    0.00533002,
	 0.00307257, -0.00509595, 1.08168,
};

// Bradford chromatic adaptation transforms between ACES white point (D60) and sRGB white point (D65)
static const float3x3 D65_2_D60_CAT =
{
	 1.01303,    0.00610531, -0.014971,
	 0.00769823, 0.998165,   -0.00503203,
	-0.00284131, 0.00468516,  0.924507,
};

// REC 709 primaries
static const float3x3 XYZ_2_sRGB_MAT =
{
	 3.2409699419, -1.5373831776, -0.4986107603,
	-0.9692436363,  1.8759675015,  0.0415550574,
	 0.0556300797, -0.2039769589,  1.0569715142,
};

static const float3x3 XYZ_2_AP1_MAT =
{
	 1.6410233797, -0.3248032942, -0.2364246952,
	-0.6636628587,  1.6153315917,  0.0167563477,
	 0.0117218943, -0.0082844420,  0.9883948585,
};

static const float3x3 XYZ_2_AP0_MAT =
{ 
	 1.0498110175, 	0.0000000000, -0.0000974845,
    -0.4959030231,	1.3733130458,  0.0982400361,
	 0.0000000000, 	0.0000000000,  0.9912520182
};

static const float3x3 sRGB_2_XYZ_MAT =
{
	 0.4124564, 0.3575761, 0.1804375,
	 0.2126729, 0.7151522, 0.0721750,
	 0.0193339, 0.1191920, 0.9503041,
};

static const float3x3 AP0_2_AP1_MAT =
{
	 1.4514393161, -0.2365107469, -0.2149285693,
	-0.0765537734,  1.1762296998, -0.0996759264,
	 0.0083161484, -0.0060324498,  0.9977163014,
};


static const float3x3 AP1_2_XYZ_MAT = 
{
	 0.6624541811, 0.1340042065, 0.1561876870,
	 0.2722287168, 0.6740817658, 0.0536895174,
	-0.0055746495, 0.0040607335, 1.0103391003,
};

static const float3x3 AP1_2_AP0_MAT = 
{
   	 0.6954522414, 0.1406786965, 0.1638690622,
   	 0.0447945634, 0.8596711185, 0.0955343182,
  	-0.0055258826, 0.0040252103, 1.0015006723
};

static const float3 AP1_RGB2Y = 
{
     0.2722287168,
     0.6740817658,
     0.0536895174
};

static const float3x3 ODT_SAT_MAT = 
{
    0.949056, 0.0471857, 0.00375827,
    0.019056, 0.9771860, 0.00375827,
    0.019056, 0.0471857, 0.93375800
};

static const float3x3 RRT_SAT_MAT = 
{
    0.9708890, 0.0269633, 0.00214758,
    0.0108892, 0.9869630, 0.00214758,
    0.0108892, 0.0269633, 0.96214800
};

static const float3x3 sRGB_2_AP0 = 
{
    0.4397010, 0.3829780, 0.1773350,
    0.0897923, 0.8134230, 0.0967616,
    0.0175440, 0.1115440, 0.8707040
};

static const float3x3 sRGB_2_AP1 = 
{
    0.61319, 0.33951, 0.04737,
    0.07021, 0.91634, 0.01345,
    0.02062, 0.10957, 0.86961
};


// Transformations from RGB to other color representations
float rgb_2_hue(float3 rgb) 
{
  // Returns a geometric hue angle in degrees (0-360) based on RGB values.
  // For neutral colors, hue is undefined and the function will return a quiet NaN value.
  float hue;
  if (rgb[0] == rgb[1] && rgb[1] == rgb[2]) 
  {
    hue = 0; //RGB triplets where RGB are equal have an undefined hue
  } 
  else 
  {
    hue = (180 / PI) * atan2(sqrt(3) * (rgb[1] - rgb[2]), 2 * rgb[0] - rgb[1] - rgb[2]);
  }
    
  if (hue < 0.0) 
  		hue += 360.0;

  return hue;
}

float rgb_2_yc(float3 rgb, float ycRadiusWeight = 1.75)
{
  // Converts RGB to a luminance proxy, here called YC
  // YC is ~ Y + K * Chroma
  // Constant YC is a cone-shaped surface in RGB space, with the tip on the 
  // neutral axis, towards white.
  // YC is normalized: RGB 1 1 1 maps to YC = 1
  //
  // ycRadiusWeight defaults to 1.75, although can be overridden in function 
  // call to rgb_2_yc
  // ycRadiusWeight = 1 -> YC for pure cyan, magenta, yellow == YC for neutral 
  // of same value
  // ycRadiusWeight = 2 -> YC for pure red, green, blue  == YC for  neutral of 
  // same value.

    float r = rgb[0]; 
    float g = rgb[1]; 
    float b = rgb[2];
  
    float chroma = sqrt(b * (b - g) + g * (g - r) + r * (r - b));

    return (b + g + r + ycRadiusWeight * chroma) / 3;
}

float rgb_2_saturation(float3 rgb)
{
    float minrgb = min(min(rgb.r, rgb.g), rgb.b);
    float maxrgb = max(max(rgb.r, rgb.g), rgb.b);
    return (max(maxrgb, 1e-10) - max(minrgb, 1e-10) / max(maxrgb, 1e-2));
}

// Textbook monomial to basis-function conversion matrix.
static const float3x3 M = 
{
	 0.5, -1.0, 0.5,
	-1.0,  1.0, 0.5,
	 0.5,  0.0, 0.0
};

float segmented_spline_c5_fwd(float x)
{
  
  float coefsLow[6] = { -4.0000000000, -4.0000000000, -3.1573765773, -0.4852499958, 1.8477324706, 1.8477324706 }; 
  float coefsHigh[6] = { -0.7185482425, 2.0810307172, 3.6681241237, 4.0000000000, 4.0000000000, 4.0000000000 };
  float2 minPoint = float2(0.18 * exp2(-15), 0.0001);
  float2 midPoint = float2(0.18, 4.8);
  float2 maxPoint = float2(0.18 * exp2(18), 10000);
  float slopeLow = 0;
  float slopeHigh = 0;
	
  const int N_KNOTS_LOW = 4;
  const int N_KNOTS_HIGH = 4;

  // Check for negatives or zero before taking the log
  //float xCheck = x <= 0 ? xCheck = 0.00006103515 : x;
  float xCheck = x <= 0 ? exp2(-14) : x;
  float logx = log10(xCheck); 

  float logy;

  if (logx <= log10(minPoint.x))
  { 
		logy = logx * slopeLow + (log10(minPoint.y) - slopeLow * log10(minPoint.x));
  } 
  else if ((logx > log10(minPoint.x)) && (logx < log10(midPoint.x))) 
  {
    	
    	float knot_coord = (N_KNOTS_LOW-1) * (logx - log10(minPoint.x)) / (log10(midPoint.x) - log10(minPoint.x));
    	int j = knot_coord;
    	float t = knot_coord - j;

    	float3 cf = { coefsLow[ j], coefsLow[ j + 1], coefsLow[ j + 2]};
    
    	float3 monomials = { t * t, t, 1 };
    	logy = dot(monomials, mul(cf, M));

  } 
  else if ((logx >= log10(midPoint.x)) && (logx < log10(maxPoint.x))) 
  {
    	
    	float knot_coord = (N_KNOTS_HIGH-1) * (logx - log10(midPoint.x)) / (log10(maxPoint.x) - log10(midPoint.x));
    	int j = knot_coord;
    	float t = knot_coord - j;

    	float3 cf = { coefsHigh[ j], coefsHigh[ j + 1], coefsHigh[ j + 2]}; 

    	float3 monomials = { t * t, t, 1. };
    	logy = dot(monomials, mul(cf, M));

  }
  else
  {
    	logy = logx * slopeHigh + (log10(maxPoint.y) - slopeHigh * log10(maxPoint.x));
  }

  return pow(10, logy);
  
}

float segmented_spline_c9_fwd(float x) //48 nits
{    
      float coefsLow[10] = { -1.6989700043, -1.6989700043, -1.4779000000, -1.2291000000, -0.8648000000, -0.4480000000, 0.0051800000, 0.4511080334, 0.9113744414, 0.9113744414 };
      float coefsHigh[10] = { 0.5154386965, 0.8470437783, 1.1358000000, 1.3802000000, 1.5197000000, 1.5985000000, 1.6467000000, 1.6746091357, 1.6878733390, 1.6878733390 };
      float2 minPoint = float2(segmented_spline_c5_fwd(0.18 * exp2(-6.5)), 0.02);
      float2 midPoint = float2(segmented_spline_c5_fwd(0.18), 4.8);
      float2 maxPoint = float2(segmented_spline_c5_fwd(0.18 * exp2(6.5)), 48);
      float slopeLow = 0;
      float slopeHigh = 0.04;

  const int N_KNOTS_LOW = 8;
  const int N_KNOTS_HIGH = 8;

  // Check for negatives or zero before taking the log. 
  float xCheck = x <= 0 ? 1e-4 : x;

  float logx = log10(xCheck);
  float logy;

  if (logx <= log10(minPoint.x)) 
  { 

    	logy = logx * slopeLow + (log10(minPoint.y) - slopeLow * log10(minPoint.x));

  } 
  else if ((logx > log10(minPoint.x)) && (logx < log10(midPoint.x)))
  {

    	float knot_coord = (N_KNOTS_LOW - 1) * (logx - log10(minPoint.x)) / (log10(midPoint.x) - log10(minPoint.x));
    	int j = knot_coord;
    	float t = knot_coord - j;

    	float3 cf = { coefsLow[ j], coefsLow[ j + 1], coefsLow[ j + 2]};
    
    	float3 monomials = { t * t, t, 1 };
    	logy = dot(monomials, mul(cf, M));

  } 
  else if ((logx >= log10(midPoint.x)) && (logx < log10(maxPoint.x)))
  {

    	float knot_coord = (N_KNOTS_HIGH - 1) * (logx - log10(midPoint.x)) / (log10(maxPoint.x) - log10(midPoint.x));
   		int j = knot_coord;
    	float t = knot_coord - j;

    	float3 cf = { coefsHigh[ j], coefsHigh[ j + 1], coefsHigh[ j + 2]}; 

    	float3 monomials = { t * t, t, 1 };
    	logy = dot(monomials, mul(cf, M));

  } 
  else
  {

    	logy = logx * slopeHigh + (log10(maxPoint.y) - slopeHigh * log10(maxPoint.x));

  }

  return pow(10, logy);
  
}

// ------- Red modifier functions
float cubic_basis_shaper(float x, float w)
{
  float M[4][4] = 
  { 
  		{ -1/6,  3/6, -3/6,  1/6 },
        {  3/6, -6/6,  3/6,  0/6 },
        { -3/6,  0/6,  3/6,  0/6 },
        {  1/6,  4/6,  1/6,  0/6 } 
  };
  
  float knots[5] = { -0.5 * w, -0.25 * w, 0, 0.25 * w, 0.5 * w };
  
  float y = 0;
  if ((x > knots[0]) && (x < knots[4]))
  {  
    float knot_coord = (x - knots[0]) * 4.0 / w;  
    int j = knot_coord;
    float t = knot_coord - j;
      
    float monomials[4] = { t * t * t, t * t, t, 1.0 };

    // (if/else structure required for compatibility with CTL < v1.5.)
    if (j == 3) {
      y = monomials[0] * M[0][0] + monomials[1] * M[1][0] + 
          monomials[2] * M[2][0] + monomials[3] * M[3][0];
    } else if (j == 2) {
      y = monomials[0] * M[0][1] + monomials[1] * M[1][1] + 
          monomials[2] * M[2][1] + monomials[3] * M[3][1];
    } else if (j == 1) {
      y = monomials[0] * M[0][2] + monomials[1] * M[1][2] + 
          monomials[2] * M[2][2] + monomials[3] * M[3][2];
    } else if (j == 0) {
      y = monomials[0] * M[0][3] + monomials[1] * M[1][3] + 
          monomials[2] * M[2][3] + monomials[3] * M[3][3];
    } else {
      y = 0.0;
    }
  }
  return y * 1.5;
}

float center_hue(float hue, float centerH)
{
  	float hueCentered = hue - centerH;
  	if (hueCentered < -180)
  		hueCentered += 360;
  	else if (hueCentered > 180) 
  		hueCentered -= 360;
  	return hueCentered;
}

float sigmoid_shaper(float x)
{
    // Sigmoid function in the range 0 to 1 spanning -2 to +2.
    float t = max(1 - abs(0.5 * x), 0);
    float y = 1 + sign(x) * (1 - t * t);
    return 0.5 * y;
}

// ------- Glow module functions
float glow_fwd(float ycIn, float glowGainIn, float glowMid)
{
   float glowGainOut;

   if (ycIn <= 2 / 3 * glowMid) 
   {
    	glowGainOut = glowGainIn;
   } 
   else if (ycIn >= 2 * glowMid) 
   {
    	glowGainOut = 0;
   } 
   else 
   {
    	glowGainOut = glowGainIn * (glowMid / ycIn - 0.5);
   }

   return glowGainOut;
}

static const float RRT_GLOW_GAIN = 0.05;
static const float RRT_GLOW_MID = 0.08;

static const float RRT_RED_SCALE = 0.82;
static const float RRT_RED_PIVOT = 0.03;
static const float RRT_RED_HUE = 0.0;
static const float RRT_RED_WIDTH = 135.0;
static const float RRT_SAT_FACTOR = 0.96;

float3 RRT(float3 aces)
{
  	// --- Glow module --- //
    float saturation = rgb_2_saturation(aces);
    float ycIn = rgb_2_yc(aces);
    float s = sigmoid_shaper((saturation - 0.4) / 0.2);
    float addedGlow = 1 + glow_fwd(ycIn, RRT_GLOW_GAIN * s, RRT_GLOW_MID);
    aces *= addedGlow;

	// --- Red modifier --- //
    float hue = rgb_2_hue(aces);
    float centeredHue = center_hue(hue, RRT_RED_HUE);
    float hueWeight = cubic_basis_shaper(centeredHue, RRT_RED_WIDTH);

    aces.r += hueWeight * saturation * (RRT_RED_PIVOT - aces.r) * (1 - RRT_RED_SCALE);

    // --- ACES to RGB rendering space --- //
    aces = clamp(aces, 0, 65535); // avoids saturated colors from becoming positive in the matrix
    float3 rgbPre = mul(AP0_2_AP1_MAT, aces);                   
    rgbPre = clamp(rgbPre, 0, 65535);
    
    // --- Global desaturation --- //
    rgbPre = lerp(dot(rgbPre, AP1_RGB2Y), rgbPre, RRT_SAT_FACTOR);

    // Apply the tonescale independently in rendering-space RGB
    float3 rgbPost;
    rgbPost[0] = segmented_spline_c5_fwd(rgbPre[0]);
    rgbPost[1] = segmented_spline_c5_fwd(rgbPre[1]);
    rgbPost[2] = segmented_spline_c5_fwd(rgbPre[2]);

    // // OCES to RGB rendering space
    float3 rgbOCES = mul(AP1_2_AP0_MAT, rgbPost);
    return rgbOCES;
}

// Transformations between CIE XYZ tristimulus values and CIE x,y 
// chromaticity coordinates
float3 XYZ_2_xyY(float3 XYZ)
{  
  float3 xyY;
  float divisor = (XYZ[0] + XYZ[1] + XYZ[2]);
  if (divisor == 0) divisor = 1e-10;
  xyY[0] = XYZ[0] / divisor;
  xyY[1] = XYZ[1] / divisor;  
  xyY[2] = XYZ[1];
  
  return xyY;
}

float3 xyY_2_XYZ(float3 xyY)
{
  float3 XYZ;
  XYZ[0] = xyY[0] * xyY[2] / max(xyY[1], 1e-10);
  XYZ[1] = xyY[2];  
  XYZ[2] = (1.0 - xyY[0] - xyY[1]) * xyY[2] / max(xyY[1], 1e-10);

  return XYZ;
}

float Y_2_linCV(float Y, float Ymax, float Ymin)
{
    return (Y - Ymin) / (Ymax - Ymin);
}

//Gamma compensation for surrounding viewing area
static const float DIM_SURROUND_GAMMA = 0.9811;
float3  darkSurround_to_dimSurround(float3 linearCV)
{
    float3 XYZ = mul(AP1_2_XYZ_MAT, linearCV);                    

    float3 xyY = XYZ_2_xyY(XYZ);
    xyY[2] = clamp(xyY[2], 0, 65535); 
    xyY[2] = pow(xyY[2], DIM_SURROUND_GAMMA);
    XYZ = xyY_2_XYZ(xyY);

    return mul(XYZ_2_AP1_MAT, XYZ);                               
}

// NOTE: The EOTF is *NOT* gamma 2.4, it follows IEC 61966-2-1:1999
static const float DISPGAMMA = 2.4; 
static const float OFFSET = 0.055;

float3 ODT_2_sRGB_D65(float3 oces)
{
    float3 rgbPre = mul(AP0_2_AP1_MAT, oces);                     

    // Apply the tonescale independently in rendering-space RGB
    float3 rgbPost;
    rgbPost[0] = segmented_spline_c9_fwd(rgbPre[0]);
    rgbPost[1] = segmented_spline_c9_fwd(rgbPre[1]);
    rgbPost[2] = segmented_spline_c9_fwd(rgbPre[2]);

    const float CINEMA_WHITE = 48.0;
    const float CINEMA_BLACK = CINEMA_WHITE / 2400;
    // Scale luminance to linear code value
    float3 linearCV;
    linearCV[0] = Y_2_linCV(rgbPost[0], CINEMA_WHITE, CINEMA_BLACK);
    linearCV[1] = Y_2_linCV(rgbPost[1], CINEMA_WHITE, CINEMA_BLACK);
    linearCV[2] = Y_2_linCV(rgbPost[2], CINEMA_WHITE, CINEMA_BLACK);    

    // Apply gamma adjustment to compensate for dim surround
    linearCV = darkSurround_to_dimSurround(linearCV);

    // Apply desaturation to compensate for luminance difference
    const float ODT_SAT_FACTOR = 0.93;
    linearCV = lerp(dot(linearCV, AP1_RGB2Y), linearCV, ODT_SAT_FACTOR);
    
    // Convert to display primary encoding
    // Rendering space RGB to XYZ
    float3 XYZ = mul(AP1_2_XYZ_MAT, linearCV);    

    // Apply CAT from ACES white point to assumed observer adapted white point
    XYZ = mul(D60_2_D65_CAT, XYZ);

    // CIE XYZ to display primaries
    linearCV = mul(XYZ_2_sRGB_MAT, XYZ);

    // Handle out-of-gamut values
    // Clip values < 0 or > 1 (i.e. projecting outside the display primaries)
    linearCV = saturate(linearCV);

    return linearCV;
   
}

//	ACES sRGB D65 Output Transform
//  Input is scene-referred linear values in the sRGB gamut
//  Output is output-referred linear values in the sRGB gamut
float3 ACESOutputTransformsRGBD65(float3 c)
{
	const float3x3 sRGB_2_AP0 = mul(XYZ_2_AP0_MAT, mul(D65_2_D60_CAT, sRGB_2_XYZ_MAT));
	float3 aces = mul(sRGB_2_AP0, c);
	float3 oces = RRT(aces);
  	return ODT_2_sRGB_D65(oces);
}

float3 toneACESRef(float3 c)
{
	return c = ACESOutputTransformsRGBD65(c);
}
