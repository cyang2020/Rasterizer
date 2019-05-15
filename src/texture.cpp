#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>


namespace CGL {
    Color Texture::sample(const SampleParams &sp) {
        // Parts 5 and 6: Fill this in.
        int level = get_level(sp);
        if (sp.psm == 0) { //pnearest
            if (sp.lsm == 0){
                return sample_nearest(sp.p_uv, 0);
            }else if (sp.lsm == 1){
                return sample_nearest(sp.p_uv, level);
            }else if (sp.lsm == 2){
                int lowLevel = (int) floor(level);
                Color u0 = (level-lowLevel)*sample_nearest(sp.p_uv, lowLevel) + (lowLevel+1-level)*sample_nearest(sp.p_uv, lowLevel+1);
                return u0;
            }
        } else {
            if (sp.lsm == 0){
                return sample_bilinear(sp.p_uv, 0);
            }else if (sp.lsm == 1){
                return sample_bilinear(sp.p_uv, level);
            }else if (sp.lsm == 2) {
                int lowLevel = (int) floor(level);
                Color u0 = (level-lowLevel)*sample_bilinear(sp.p_uv, lowLevel) + (lowLevel+1-level)*sample_bilinear(sp.p_uv,lowLevel+1);
                return u0;
            }
        }
    }
    
    float Texture::get_level(const SampleParams &sp) {
        // Optional helper function for Parts 5 and 6
        Vector2D uv = Vector2D(sp.p_uv[0]*width, sp.p_uv[1]*height);
        Vector2D du = Vector2D(sp.p_dx_uv[0]*width, sp.p_dx_uv[1]*height);
        Vector2D dv = Vector2D(sp.p_dy_uv[0]*width, sp.p_dy_uv[1]*height);
        
        float L = max((du-uv).norm(),(dv-uv).norm());
        
        L = log2f(L);
        if (L<0){
            L=0;
        }else if (L>=kMaxMipLevels){
            L=kMaxMipLevels;
        }
        return L;
    }
    
    // Returns the nearest sample given a particular level and set of uv coords
Color Texture::sample_nearest(Vector2D uv, int level) {
    int width = mipmap[level].width;
    int height = mipmap[level].height;
    int u = uv[0]*width;
    int v = uv[1]*height;
    if (u>=width){
        u=width-1;
    } else if (v>=height){
        v=height-1;
    }
    int i = floor(3*(width*v+u));
    return Color(&mipmap[level].texels[i]);
}
    

    // Returns the bilinear sample given a particular level and set of uv coords
    Color Texture::sample_bilinear(Vector2D uv, int level) {
        // Optional helper function for Parts 5 and 6
        // Feel free to ignore or create your own
        int height = mipmap[level].height;
        int width = mipmap[level].width;
        
        float u = uv[0]*width;
        float v = uv[1]*height;
        if (u>=width){
            u=width-1;
        } else if (v>=height){
            v=height-1;
        }
        int floorU = floor(u);
        int topU = floorU+1;
        int floorV = floor(v);
        int topV = floorV+1;
        float du = u + (-1) * floorU;
        float dv = v + (-1) * floorV;
        if (topU >= width || topV >= height){
            if (topU >= width && topV >= height){
                return (Color(&mipmap[level].texels[floor(3*(floorU+floorV*width))]));
            }else if (topU >= width){
                Color u00 = Color(&mipmap[level].texels[floor(3*(floorU+floorU*width))]);
                Color u01 = Color(&mipmap[level].texels[floor(3*(floorU+topV*width))]);
                Color u0 = u00 + dv*(u01 + -1 * u00);
                return u0;
            }else{
                Color u00 = Color(&mipmap[level].texels[floor(3*(floorU+floorV*width))]);
                Color u10 = Color(&mipmap[level].texels[floor(3*(topU+floorV*width))]);
                Color u0 = u00 + du*(u10 + -1 * u00);
                return u0;
            }
        }else{
            
            Color u00 = Color(&mipmap[level].texels[floor(3*(floorU+floorV*width))]);
            Color u01 = Color(&mipmap[level].texels[floor(3*(floorU+topV*width))]);
            Color u10 = Color(&mipmap[level].texels[floor(3*(topU+floorV*width))]);
            Color u11 = Color(&mipmap[level].texels[floor(3*(topU+topV*width))]);
            
            Color u0 = u00 + du*(u10 + -1 * u00);
            Color u1 = u01 + du*(u11 + -1 * u01);
            Color uf = u0 + dv*(u1 + -1 * u0);
            return uf;
        }
    }


/****************************************************************************/



inline void uint8_to_float(float dst[3], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[3]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
}

void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    //assert (width > 0);
    height = max(1, height / 2);
    //assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(3 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
    int currLevelPitch = currLevel.width * 3; // 32 bit RGB

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[3];
    float input[3];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      //assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (3 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      //assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        3 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
        }
      }
    }
  }
}

}
