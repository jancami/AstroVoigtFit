from functions.voigtprofile import voigt_profile

def test_voigt():
  assert voigt_profile(2.5,3,9) == 0.030695779845761356
