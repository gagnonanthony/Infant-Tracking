FROM scilus/scilus-flows:1.5.0

LABEL version="ChildBrainFlow-1.0.0"

RUN wget https://www.dropbox.com/scl/fi/izy3yie5cdqqi7my2ogfa/FS_BN_GL_SF_utils.tar.gz?rlkey=kt8qi4ye7h3om3ymnq0n0h570&dl=0 -O FS_BN_GL_SF_utils.tar.gz && \
    tar -xvf FS_BN_GL_SF_utils.tar && \
    rm FS_BN_GL_SF_utils.tar

# Installing freesurfer on top of scilus:1.5.0
WORKDIR /root
RUN wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-centos7_x86_64-7.2.0.tar.gz -O fs.tar.gz && \
    tar --no-same-owner -xzvf fs.tar.gz && \
    mv freesurfer /usr/local && \
    rm fs.tar.gz
COPY license.txt /usr/local/freesurfer/license.txt
RUN apt-get update && \
    apt-get install csh tcsh
RUN wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2 \
    tar xjf parallel-latest.tar.bz2 \
    cd parallel-* && ./configure && make && make install \
    echo 'will cite' | parallel --citation 1> /dev/null 2> /dev/null &

# Setup freesurfer env
ENV OS=Linux
ENV PATH=${PATH}:/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV FREESURFER_HOME=/usr/local/freesurfer
ENV FREESURFER=/usr/local/freesurfer
ENV SUBJECTS_DIR=/usr/local/freesurfer/subjects
ENV LOCAL_DIR=/usr/local/freesurfer/local
ENV FSFAST_HOME=/usr/local/freesurfer/fsfast
ENV FMRI_ANALYSIS_DIR=/usr/local/freesurfer/fsfast
ENV FUNCTIONALS_DIR=/usr/local/freesurfer/sessions
ENV FS_LICENSE=/usr/local/freesurfer/license.txt

# set default fs options
ENV FS_OVERRIDE=0
ENV FIX_VERTEX_AREA=""
ENV FSF_OUTPUT_FORMAT=nii.gz

# mni env requirements
ENV MINC_BIN_DIR=/usr/local/freesurfer/mni/bin
ENV MINC_LIB_DIR=/usr/local/freesurfer/mni/lib
ENV MNI_DIR=/usr/local/freesurfer/mni
ENV MNI_DATAPATH=/usr/local/freesurfer/mni/data
ENV MNI_PERL5LIB=/usr/local/freesurfer/mni/share/perl5
ENV PERL5LIB=/usr/local/freesurfer/mni/share/perl5