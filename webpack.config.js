const webpack = require('webpack');

module.exports = {
    //插件项
    //页面入口文件配置
    entry: {
        index: './public/js/index.js'
    },
    //入口文件输出配置
    output: {
        path: 'public/javascripts',
        filename: '[name].js'
    },
    module: {
        //加载器配置
        loaders: [ {
            test: /\.css$/,
            loader: 'style-loader!css-loader'
        }, {
            test: /\.scss$/,
            loader: 'style!css!sass?sourceMap'
        }, {
            test: /\.(png|jpg)$/,
            loader: 'url-loader?limit=8192'
        }]
    },
    //其它解决方案配置
    resolve: {
        extensions: ['', '.js', '.json', '.scss'],
        alias: {}
    }
};