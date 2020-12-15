#!groovy
pipeline {
    agent any

    stages {
        stage('Build') {
            agent{
                image 'braintwister/ubuntu-18.04-gcc-6:0.3'
            }
            steps {
            sh 'cmake -Bbuild -H. && cd build && make'
            }
        }
        stage('Test') {
            steps {
                sh 'cd build && make test'
            }
        }
    }
}
