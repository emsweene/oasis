library(hexSticker)
library(desc)
desc = desc::description$new()
package = desc$get("Package")
outline = "#e9371c"
background = "lightgray"
sticker("icon.png",	
        package = package,
        h_fill = background,
        h_color = outline, 
        s_width = 0.45, 
        s_height = 0.9,
        s_x = 1,
        filename = "sticker.png")

usethis::use_build_ignore(
  c("icon.png", "sticker.R", "sticker.png"))

