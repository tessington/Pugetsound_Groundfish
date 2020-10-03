
install.packages(c("rnaturalearth", "rnaturalearthdata"))
                 

require(ggplot2)
require(cowplot)
require(sf)        
require(rnaturalearth)
require(rnaturalearthdata)
require(rnaturualearthhires)

theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-140.15, -100.00), ylim = c(20.00, 55.00), expand = FALSE) +
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size = 18)) +
  scale_x_continuous(breaks = c(-140, -120, -100)) +
  scale_y_continuous(breaks = c(25, 35, 45, 55))


        
        
        
ggsave("wcoast.pdf", height = 5, width = 5)

