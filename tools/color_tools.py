import pickle, os, sys

# transparency:
transparency2hex = {}
transparency2hex[100] = 'FF'
transparency2hex[95] = 'F2'
transparency2hex[90] = 'E6'
transparency2hex[85] = 'D9'
transparency2hex[80] = 'CC'
transparency2hex[75] = 'BF'
transparency2hex[70] = 'B3'
transparency2hex[65] = 'A6'
transparency2hex[60] = '99'
transparency2hex[55] = '8C'
transparency2hex[50] = '80'
transparency2hex[45] = '73'
transparency2hex[40] = '66'
transparency2hex[35] = '59'
transparency2hex[30] = '4D'
transparency2hex[25] = '40'
transparency2hex[20] = '33'
transparency2hex[15] = '26'
transparency2hex[10] = '1A'
transparency2hex[5] = '0D'
transparency2hex[0] = '00'

# get path of this script, assume tools are in directory 'tools/' in the parent directory of this code
type2size2name2rgblist = pickle.load(open(os.path.abspath(os.path.dirname(__file__))+'/support_files/brewer-palettes.pickled', 'rb'))

def transparency_to_hex(T):
	return transparency2hex[T]
	
	
###some methods to get color gradients
#  --------- SOURCE http://bsou.io/posts/color-gradients-with-python ------------
def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  # Pass 16 to the integer function for change of base
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]


def RGB_to_hex(RGB,T=0):
  ''' [255,255,255] -> "#FFFFFF" '''
  # Components need to be integers for hex to make sense
  RGB = [int(x) for x in RGB]
  if T == 0:
    return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])
  else:
    thex = transparency2hex[T]
    return "#"+thex+"".join(["0{0:x}".format(v) if v < 16 else
              "{0:x}".format(v) for v in RGB])


def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return  color_dict(RGB_list)['hex']


  def linear_gradient_3colors(start_hex, middle_hex, finish_hex="#FFFFFF", n=10):
    ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''

    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    m = hex_to_RGB(middle_hex)
    f = hex_to_RGB(finish_hex)

    middle_point = int(n/2.0)
    if n%2 != 0: print('Warning number of colors is not even, will take',middle_point,'as the middle point')

    #from start to middle:
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, middle_point):
      # Interpolate RGB vector for color at the current value of t
      curr_vector = [
        int(s[j] + (float(t)/(n-1))*(m[j]-s[j]))
        for j in range(3)
      ]
      # Add it to our list of output colors
      RGB_list.append(curr_vector)

    # from middle to finish
    for t in range(middle_point, n):
      # Interpolate RGB vector for color at the current value of t
      curr_vector = [
        int(m[j] + (float(t)/(n-1))*(f[j]-m[j]))
        for j in range(3)
      ]
      # Add it to our list of output colors
      RGB_list.append(curr_vector)

    return  color_dict(RGB_list)['hex']


def get_brewer_colorlist_in_hex(ptype, N, transparency, name = None):
  colorlist_rgb = []

  if ptype == 'qual':
    if name == None:
      if transparency > 10:
        name = 'dark2'
      else:
        name = 'set1'

    colorlist_rgb = type2size2name2rgblist[ptype][str(N)][name]

  else:
    if name == None:
      if ptype == 'div':
        name = 'spectral'
      else:
        name = 'ylorrd'

    if N <= 7:
      colorlist_rgb = type2size2name2rgblist[ptype][str(N+2)][name][2:]
    elif N == 8:
      colorlist_rgb = type2size2name2rgblist[ptype]['9'][name][1:]
    else:
      colorlist_rgb = type2size2name2rgblist[ptype]['9'][name]

  colorlist = [RGB_to_hex(RGB, T = transparency) for RGB in colorlist_rgb]
  return colorlist



