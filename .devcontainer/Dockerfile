FROM condaforge/miniforge3:latest

ARG USERNAME=pycharm
ARG USER_UID=1000
ARG USER_GID=1000

# System packages + Node.js
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    zsh \
    sudo \
    nodejs \
    npm \
    && rm -rf /var/lib/apt/lists/*

# Rename existing ubuntu user (UID 1000) to pycharm
RUN usermod -l $USERNAME -d /home/$USERNAME -m ubuntu \
    && groupmod -n $USERNAME ubuntu \
    && chsh -s /bin/zsh $USERNAME \
    && echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/$USERNAME

# Give user ownership of conda
RUN chown -R $USERNAME:$USERNAME /opt/conda

USER $USERNAME
WORKDIR /home/$USERNAME

# Initialize conda for the user's shells
RUN conda init bash && conda init zsh

WORKDIR /workspaces/project