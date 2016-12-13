const webpack = require('webpack')
const ExtractTextPlugin = require('extract-text-webpack-plugin')
require('sass-loader')
module.exports = {
    //插件项
    //页面入口文件配置
    entry: {
        index: './public/js/index.js',
        article:'./public/js/article.js',
        
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
            loader: ExtractTextPlugin.extract("style", 'css!sass')
        }, {
            test: /\.(png|jpg)$/,
            loader: 'url-loader?limit=8192'
        }]
    },
    plugins: [
        new ExtractTextPlugin("/public/stylesheets/[name].css") 
    ],
    //其它解决方案配置
    resolve: {
        extensions: ['', '.js', '.json', '.scss'],
        alias: {}
    }
    
};