
.get_colors_25 = function() {
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  return(c25)
}

.get_colors_12 = function() {
  c12 =  c(
    "dodgerblue2", "#E31A1C",
    "green4","deeppink1",
    "#6A3D9A", "#FF7F00",
    "black", "gold1",
    "palegreen2","darkorange4",
    "orchid1", "darkturquoise")
  return(c12)
}

.get_tirtl_pallette = function() {
  TIRTL_pallette = c("#1D5F8A","#E76F47","#A53828","#531E1E","#116D4D","#ECC30B","#DA6140","#BABAA0","#38BBB7","#509E6E","grey90")
  return(TIRTL_pallette)
}

.tirtl_colors_gradient = function(
    palette = c("viridis","sea_green","sea1", "sea_multicolor"),
    n=20,
    verbose = FALSE
    ) {
  default = "viridis"
  if(is.null(palette)) palette = default
  palette = palette[1]
  palettes = c("viridis","sea_green","sea1", "sea_multicolor")
  if(!palette %in% palettes) palette = default

  ### Color gradients

  viridis = function(n) {
    grDevices::hcl.colors(palette = "Viridis",n=n, rev = TRUE)
  }

  sea1 <- c(
    "#001f23",  # Deep sea teal (very dark)
    "#014d64",  # Ocean Deep
    "#2e8c85",  # Tidepool Teal
    "#49a9a1",  # Sea Glass
    "#79d4c4",  # Seafoam
    "#b0ebe8"   # Aqua Light (brightest tint)
  ) %>% rev()

  sea_green <- c(
    "#003220",  # Deep algae green (very dark)
    "#1b5e4a",  # Turtle green
    "#2e8c85",  # Seaweed teal
    "#49a97a",  # Marine green
    "#79d4a2",  # Light green-seafoam
    "#b0f0c2"   # Pale seafoam green
  ) %>% rev()

  sea_multicolor <- #colorRampPalette(
    c(
    "#001219",  # Abyssal Ink (deep ocean)
    "#014d64",  # Ocean Deep (dark teal)
    "#2e8c85",  # Tidepool Teal
    "#79d4c4",  # Seafoam
    "#b0ebe8",  # Aqua Light
    "#6a8e3a",  # Seagrass Green
    "#a2c523",  # Sunlit Algae
    "#f2c28b",  # Coral Sand
    "#f4b6b6",  # Shell Pink
    "#b49d79",  # Driftwood
    "#5b4e2d"   # Turtle Shell (earthy end)
  ) %>% rev()

  #col_fun = get(palette)
  #cols = do.call(col_fun, args = list(n))
  if(palette == "viridis") {
    cols = viridis(n)
  } else {
    cols = get(palette)
    if(n>length(cols)) {
      if(verbose) warning(paste("Returning all", length(cols), "colors"))
      n=length(cols)
    }
    cols = cols[1:n]
  }
  return(cols)
}

.tirtl_colors_distinct = function(
    palette = c("tirtl","sea", "sea_alt", "turtle"),
    n=Inf,
    verbose = FALSE
    ) {
  default = "tirtl"  ## use sea as default
  if(is.null(palette)) palette = default
  palette = palette[1]
  palettes = c("tirtl", "sea", "sea_alt", "turtle")
  if(!palette %in% palettes) palette = default

  tirtl = c("#1D5F8A","#E76F47","#A53828","#531E1E","#116D4D","#ECC30B","#DA6140","#BABAA0","#38BBB7","#509E6E","grey90")

  ### Distinct color palettes
  sea <- c(
    "#5b4e2d",  # Turtle Shell (earthy, dark brown)
    "#49a9d9",  # Lagoon Blue (bright blue)
    "#6a8e3a",  # Seagrass Green (earthy green)
    "#fa7268",  # Reef Coral (vivid red-orange)
    "#b0ebe8",  # Aqua Light (bright cyan)
    "#395c3b",  # Kelp Forest (deep green)
    "#f2c28b",  # Coral Sand (warm beach tone)
    "#2e8c85",  # Tidepool Teal (blue-green)
    "#b24030",  # Crustacean Red (burnt red)
    "#79d4c4",  # Seafoam (pale green-cyan)
    "#78804b",  # Olive Turtle (neutral olive)
    "#014d64",  # Ocean Deep (dark teal)
    "#f4b6b6",  # Shell Pink (soft pink)
    "#a2c523",  # Sunlit Algae (bright yellow-green)
    "#cfd8d7",  # Foam Grey (cool light gray)
    "#b49d79",  # Driftwood (desaturated beige)
    "#1c5d99",  # Marine Blue (rich mid blue)
    "#e6d2ae",  # Soft Sand (light tan)
    "#001219",  # Abyssal Ink (near-black blue)
    "#f5f3e7"   # Shell White (warm off-white)
  )

  sea_alt = c(
    "#014d64",  # Ocean Deep
    "#79d4c4",  # Seafoam
    "#5b4e2d",  # Turtle Shell
    "#f2c28b",  # Coral Sand
    "#6a8e3a",  # Seagrass Green
    "#fa7268",  # Reef Coral
    "#49a9d9",  # Lagoon Blue
    "#b49d79",  # Driftwood
    "#f5f3e7",  # Shell White
    "#2e8c85",  # Tidepool Teal
    "#78804b",  # Olive Turtle
    "#1c5d99",  # Marine Blue
    "#395c3b",  # Kelp Forest
    "#e6d2ae",  # Soft Sand
    "#b0ebe8",  # Aqua Light
    "#f4b6b6",  # Shell Pink
    "#001219",  # Abyssal Ink
    "#b24030",  # Crustacean Red
    "#cfd8d7",  # Foam Grey
    "#a2c523"   # Sunlit Algae
  )

  turtle = c(
    "#beab5c",
    "#87c529",
    "#b4ccd5",
    "#a0ad61",
    "#25d8b9",
    "#6c4f3a",
    "#346c04",
    "#e9dfae",
    "#222d4d",
    "#e2600f",
    "#4e5930",
    "#5d5b6f",
    "#64bc51",
    "#94d092",
    "#9ba91a",
    "#7e7e77",
    "#8c9cbc"
  )

  trek1 = c("#8B799CFF", "#3C999CFF", "#C86C18FF", "#E2ED50FF",
            "#13A4EBFF", "#FFF7A3FF", "#944D40FF", "#524559FF",
            "#CA480DFF", "#2F7270FF", "#F9AB3CFF",
            "#2E7BC5FF", "#BFCAFEFF", "#66FFFFFF", "#B46356FF",
            "#7A4B42FF", "#D78017FF", "#8BEAFFFF",
            "#9B5928FF", "#A1B3E2FF", "#FFE705FF")

  cols25 = c("#1F78C8", "#ff0000", "#33a02c", "#6A33C2", "#ff7f00", "#565656",
             "#FFD700", "#a6cee3", "#FB6496", "#b2df8a", "#CAB2D6", "#FDBF6F",
             "#999999", "#EEE685", "#C8308C", "#FF83FA", "#C814FA", "#0000FF",
             "#36648B", "#00E2E5", "#00FF00", "#778B00", "#BEBE00", "#8B3B00",
             "#A52A3C")

  cols = get(palette)
  if(n>length(cols)) {
    if(verbose) warning(paste("Returning all", length(cols), "colors"))
    n=length(cols)
  }
  return(cols[1:n])
}
