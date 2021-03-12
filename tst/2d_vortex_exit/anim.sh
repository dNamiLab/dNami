#ffmpeg -r 8 -i ./pics/rho%04d.png -c:v libx264 -r 8 -pix_fmt yuv420p top_right.mp4 
echo 'Enter filename:'
read fname
ffmpeg -r 8 -i ./pics/rho%04d.png -c:v libx264 -r 8 -pix_fmt yuv420p $fname 
