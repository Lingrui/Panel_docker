FROM 192.168.6.28:5000/weixuan/ubuntubase:v0.2
MAINTAINER cailingrui,cailingrui@micro-helix.com
RUN apt-get install -y ImageMagick
ADD r.install.sh /
RUN sh /r.install.sh
ADD docker_panel.pl /
ADD ./bin/* /weixuan/panel/bin/
ADD ./annovar/* /weixuan/panel/annovar/
#RUN useradd mh01
WORKDIR /
CMD perl /docker_panel.pl
