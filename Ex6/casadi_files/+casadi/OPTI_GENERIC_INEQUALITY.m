function v = OPTI_GENERIC_INEQUALITY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = casadiMEX(0, 131);
  end
  v = vInitialized;
end
